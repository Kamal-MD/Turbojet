"""
Read the exported STEP files back and check them against the model.

Called by combustor_model.py after export (result in output/verify.txt), or
run on its own:  python verify_step.py

1. 11_main_liner.STEP: every drilled hole is found from its cylindrical face;
   rows are grouped by (x, d) and compared with the hole table (count, d, x,
   first angle). Proves the STEP, not just the in-memory model, is right.
2. assembly.STEP: every pair of solids is intersected. Any overlap with
   volume is reported; touching pairs (gasket faces, welded joints, bolts in
   their nuts) are listed as contacts.
"""

import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import cadquery as cq

TOL = 1e-3


def read_liner_holes(path):
    """{(x, d): [angles]} for every radial hole in the liner STEP."""
    liner = cq.importers.importStep(str(path)).val()
    holes = set()
    for f in liner.Faces():
        if f.geomType() != "CYLINDER":
            continue
        cyl = f._geomAdaptor().Cylinder()
        if abs(cyl.Axis().Direction().X()) > 0.5:
            continue   # liner bore / OD run along X: not a drilled hole
        c = f.Center()   # lies on the hole axis, inside the wall
        theta = round(math.degrees(math.atan2(c.y, c.z)) % 360, 4) % 360
        holes.add((round(c.x, 3), round(2 * cyl.Radius(), 3), theta))
    rows = defaultdict(list)
    for x, d, theta in holes:
        rows[(x, d)].append(theta)
    return rows


def name_solids(solids, parts):
    """Match imported solids to part names by volume and centre (the STEP
    round trip moves both slightly, hence the loose tolerances)."""
    ref = [(n, s) for n, wp in parts.items() for s in wp.solids().vals()]
    names = []
    for s in solids:
        hit = [n for n, r in ref
               if abs(r.Volume() / s.Volume() - 1) < 1e-3 and (r.Center() - s.Center()).Length < 0.05]
        names.append(hit[0] if hit else "?")
    return names


def verify_outputs(out_dir, rows, bolt_angles, bolt_hole_d, x_bolt, parts):
    out_dir = Path(out_dir)
    lines = ["STEP read-back: 11_main_liner.STEP",
             f"{'row':11} {'count':>5} {'d':>6} {'x':>7} {'first angle':>11}   expected (count, d, x, angle)"]
    expected = [(r, n, d, x, ph % 360) for (r, _, n, d, x, ph) in rows]
    expected.append(("bolt", len(bolt_angles), bolt_hole_d, x_bolt, min(a % 360 for a in bolt_angles)))
    found = read_liner_holes(out_dir / "11_main_liner.STEP")
    all_ok = True
    for (r, n, d, x, ph) in sorted(expected, key=lambda e: e[3]):
        key = next((k for k in found if abs(k[0] - x) < TOL and abs(k[1] - d) < TOL), None)
        if key is None:
            lines.append(f"{r:11} {'-':>5} {'-':>6} {'-':>7} {'-':>11}   {n}, {d}, {x}, {ph:.4f}  MISSING")
            all_ok = False
            continue
        ths = found.pop(key)
        first = min(ths)
        ok = len(ths) == n and abs(first - ph) < TOL
        all_ok &= ok
        lines.append(f"{r:11} {len(ths):5d} {key[1]:6.2f} {key[0]:7.2f} {first:11.4f}   "
                     f"{n}, {d}, {x}, {ph:.4f}  {'MATCH' if ok else 'MISMATCH'}")
    for (x, d), ths in found.items():
        lines.append(f"unexpected hole group x={x} d={d} count={len(ths)}")
        all_ok = False
    lines.append("RESULT: all rows match" if all_ok else "RESULT: MISMATCH")

    solids = cq.importers.importStep(str(out_dir / "assembly.STEP")).solids().vals()
    names = name_solids(solids, parts)
    overlaps, contacts = [], []
    for (i, a), (j, b) in combinations(enumerate(solids), 2):
        if a.distance(b) > 1e-6:
            continue
        v = a.intersect(b).Volume()
        (overlaps if v > 1e-4 else contacts).append((names[i], names[j], v))
    lines += ["", f"Interference: assembly.STEP, {len(solids)} solids ({names.count('?')} unmatched)"]
    if overlaps:
        lines += [f"  OVERLAP {a} / {b}: {v:.4f} mm3" for a, b, v in overlaps]
    else:
        lines.append("  no overlapping solids (volume > 1e-4 mm3)")
    lines.append("  touching pairs (contact, zero overlap):")
    lines += [f"    {a:18} - {b}" for a, b, _ in contacts]
    return "\n".join(lines)


if __name__ == "__main__":
    import combustor_model as m
    parts = {n: build() for n, build in m.PART_BUILDERS}
    print(verify_outputs(m.OUT, m.ROWS, m.BOLT_ANGLES, m.bolt_hole_d, m.x_bolt, parts))
