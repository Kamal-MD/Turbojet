"""
Parametric reference model of the propane turbojet combustor (ground test only).

Builds every part of the can combustor as a CadQuery solid in one global frame,
runs geometric / flow sanity checks and exports STEP files for use as a
reference in Fusion 360.

Frame (all dimensions in mm):
    X = engine axis, flow goes toward +X (head end -> turbine end)
    x = 0 at the inner face of the head plate (= the blind flange)
    the compressor inlet stub points along +Z
    hole angles are measured about X, from +Z toward +Y

Run:  python combustor_model.py        (outputs go to ./output/)
"""

import csv
import math
from itertools import combinations
from pathlib import Path

import cadquery as cq

OUT = Path(__file__).resolve().parent / "output"

# =============================================================================
# PARAMETERS
# -----------------------------------------------------------------------------
# Every number the model uses is registered here with unit, description and
# kind, so parameters.csv is a complete record of the build:
#   input    - value given in the design spec (do not change without review)
#   derived  - computed from inputs
#   assumed  - placeholder chosen for the model because the spec gives no
#              value; review these before ordering material
#   standard - value read from a standard (EN 1092-1 / ISO 4032 ...), verify
#              against the supplier's catalogue
#   note     - text only
# =============================================================================
PARAMS = {}


def p(name, value, unit, desc, kind="input"):
    PARAMS[name] = (value, unit, desc, kind)
    return value


# --- Wall / liner -------------------------------------------------------------
wall = p("wall", 2.5, "mm", "Wall of liner, ring and housing tube (stick-weld limit: thinner stock burns through)")
liner_OD = p("liner_OD", 76.1, "mm", "Main liner tube OD (common 76.1 x 2.5 tube)")
liner_ID = p("liner_ID", liner_OD - 2 * wall, "mm", "Main liner ID", "derived")
ring_OD = p("ring_OD", 63.5, "mm", "Head ring tube OD (sits inside the liner, forms the film slot)")
ring_ID = p("ring_ID", ring_OD - 2 * wall, "mm", "Head ring ID", "derived")
ring_L = p("ring_L", 30.0, "mm", "Head ring length from the head plate; the cooling film exits at x = ring_L")
tube_gasket_t = p("tube_gasket_t", 1.0, "mm", "Graphite gasket between liner end and head plate (placeholder)")
liner_start = p("liner_start", tube_gasket_t, "mm", "x of the liner upstream end (sits on the tube gasket)", "derived")
overlap_L = p("overlap_L", ring_L - liner_start, "mm", "Axial overlap of liner over head ring (film path length)", "derived")
film_gap = p("film_gap", (liner_ID - ring_OD) / 2, "mm", "Radial film-cooling slot between ring OD and liner ID", "derived")
liner_L = p("liner_L", 250.0, "mm", "x of the liner downstream end (liner length measured from x = 0)")
main_L = p("main_L", liner_L - liner_start, "mm", "Physical length of the main liner tube", "derived")

# --- Ring bolts (liner -> head ring) --------------------------------------------
n_bolt = p("n_bolt", 3, "-", "Number of M4 bolts holding the liner on the head ring")
bolt_d = p("bolt_d", 4.0, "mm", "Bolt nominal diameter (M4)")
bolt_hole_d = p("bolt_hole_d", 4.5, "mm", "Bolt clearance hole, match-drilled through liner + ring together")
x_bolt = p("x_bolt", 9.0, "mm", "x of the bolt row (head inside the ring, nut outside the liner)")
phase_bolt = p("phase_bolt", 0.0, "deg", "Angle of the first bolt from +Z", "assumed")
bolt_head_AF = p("bolt_head_AF", 7.0, "mm", "M4 hex head across flats (ISO 4017)", "standard")
bolt_head_k = p("bolt_head_k", 2.8, "mm", "M4 hex head height (ISO 4017)", "standard")
nut_AF = p("nut_AF", 7.0, "mm", "M4 nut across flats (ISO 4032)", "standard")
nut_m = p("nut_m", 3.2, "mm", "M4 nut height (ISO 4032)", "standard")
spacer_OD = p("spacer_OD", 10.0, "mm", "Spacer (film-gap distance piece) OD, one per bolt", "assumed")

# --- Liner holes (all radial, drilled straight, equally spaced per row) --------
n_cool = p("n_cool", 28, "-", "Cooling holes: count")
d_cool = p("d_cool", 2.0, "mm", "Cooling holes: diameter")
x_cool = p("x_cool", 18.0, "mm", "Cooling holes: x (liner wall over the ring -> feed the film slot, film exits at x = ring_L)")
phase_cool = p("phase_cool", 0.0, "deg", "Cooling holes: angle of first hole from +Z", "assumed")
n_prim = p("n_prim", 10, "-", "Primary holes: count")
d_prim = p("d_prim", 6.5, "mm", "Primary holes: diameter")
x_prim = p("x_prim", 62.0, "mm", "Primary holes: x")
phase_prim = p("phase_prim", 0.0, "deg", "Primary holes: angle of first hole from +Z", "assumed")
n_sec = p("n_sec", 14, "-", "Secondary holes: count")
d_sec = p("d_sec", 4.5, "mm", "Secondary holes: diameter")
x_sec = p("x_sec", 115.0, "mm", "Secondary holes: x")
phase_sec = p("phase_sec", 360.0 / n_sec / 2, "deg", "Secondary holes: half a secondary pitch off the primary reference row", "derived")
n_dil_row = p("n_dil_row", 7, "-", "Dilution holes per row (2 staggered rows, 14 total)")
d_dil = p("d_dil", 10.0, "mm", "Dilution holes: diameter")
x_dil1 = p("x_dil1", 145.0, "mm", "Dilution row 1: x")
x_dil2 = p("x_dil2", 170.0, "mm", "Dilution row 2: x")
phase_dil1 = p("phase_dil1", 2.5714, "deg", "Dilution row 1: chosen to maximise separation for n_prim=10, n_sec=14, n_dil_row=7; recheck if any count changes")
phase_dil2 = p("phase_dil2", phase_dil1 + 360.0 / n_dil_row / 2, "deg", "Dilution row 2: staggered half a pitch from row 1", "derived")

# --- Hole area targets (from the air split, Cd 0.6) ----------------------------
A_cool_target = p("A_cool_target", 87.1, "mm2", "Cooling hole area target")
A_prim_target = p("A_prim_target", 329.0, "mm2", "Primary hole area target")
A_sec_target = p("A_sec_target", 219.4, "mm2", "Secondary hole area target")
A_dil_target = p("A_dil_target", 1105.5, "mm2", "Dilution hole area target")
area_tol = p("area_tol", 0.03, "-", "Allowed relative deviation of drilled area from target")

# --- Housing --------------------------------------------------------------------
housing_OD = p("housing_OD", 108.0, "mm", "Housing tube OD (bolted assembly, not welded to the liner; graphite gaskets)")
housing_ID = p("housing_ID", 103.0, "mm", "Housing tube ID")

# --- Head section -----------------------------------------------------------------
head_OD = p("head_OD", 139.7, "mm", "Head pipe OD (DN125)")
head_wall = p("head_wall", 2.0, "mm", "Head pipe wall (option 3.0)")
head_ID = p("head_ID", head_OD - 2 * head_wall, "mm", "Head pipe ID", "derived")
step_t = p("step_t", 4.0, "mm", "Flat step plate between DN125 head pipe and 108 housing (spec range 3-4, NOT a cone)", "assumed")
stub_OD = p("stub_OD", 50.0, "mm", "Inlet stub OD (fits the existing 50.5 mm ID silicone hose)")
stub_ID = p("stub_ID", stub_OD - 2 * wall, "mm", "Inlet stub ID (compressor delivery)", "derived")
stub_protrude = p("stub_protrude", 40.0, "mm", "Stub height above the head pipe OD (to the compressor connection)", "assumed")

# --- Flanges: EN 1092-1 PN6 DN125 -------------------------------------------------
flange_OD = p("flange_OD", 240.0, "mm", "Flange OD (EN 1092-1 PN6 DN125)")
flange_BCD = p("flange_BCD", 200.0, "mm", "Flange bolt circle diameter")
n_flange_bolt = p("n_flange_bolt", 8, "-", "Flange bolts (M16)")
flange_hole_d = p("flange_hole_d", 18.0, "mm", "Flange bolt hole diameter")
flange_phase = p("flange_phase", 180.0 / n_flange_bolt, "deg", "Bolt holes straddle the vertical (stub) centreline", "derived")
flange_t = p("flange_t", 15.0, "mm", "Flange thickness, blind and slip-on (placeholder until the bought flange is measured)")
flange_bore = p("flange_bore", 141.5, "mm", "Slip-on (type 01) flange bore for 139.7 pipe", "standard")
fl_gasket_t = p("fl_gasket_t", 2.0, "mm", "Graphite flange gasket thickness (placeholder)")
slipon_setback = p("slipon_setback", head_wall, "mm", "Pipe end set back inside the slip-on flange face (room for the internal fillet weld)", "assumed")
hw_clear = p("hw_clear", 19.0, "mm", "Clearance behind slip-on flange for M16 nuts + spanner, up to the stub edge")
stub_x = p("stub_x", fl_gasket_t + flange_t + hw_clear + stub_OD / 2, "mm", "x of the stub axis = fl_gasket_t + flange_t + hw_clear + stub_OD/2", "derived")
head_L = p("head_L", stub_x + stub_OD / 2 + 10, "mm", "Head pipe end / step plate upstream face = stub_x + stub_OD/2 + 10", "derived")

# --- Turbine end slip joint ---------------------------------------------------------
slip_axial = p("slip_axial", 5.0, "mm", "Cold axial clearance, liner end -> end plate (liner grows ~3 mm hot)")
slip_radial = p("slip_radial", 0.5, "mm", "Radial clearance liner OD -> slip sleeve ID")
sleeve_ID = p("sleeve_ID", liner_OD + 2 * slip_radial, "mm", "Slip sleeve ID", "derived")
sleeve_OD = p("sleeve_OD", sleeve_ID + 2 * wall, "mm", "Slip sleeve OD", "derived")
sleeve_L = p("sleeve_L", 20.0, "mm", "Slip sleeve length (cold engagement = sleeve_L - slip_axial)", "assumed")
end_t = p("end_t", step_t, "mm", "End plate thickness (same plate stock as the step plate)", "assumed")
end_hole_d = p("end_hole_d", liner_ID, "mm", "End plate gas outlet bore (exit transition to turbine not modelled)", "derived")
x_end_in = p("x_end_in", liner_L + slip_axial, "mm", "x of the end plate inner face", "derived")

# --- Flow / thermal data used by the checks (design log 2026-10-14) ---------------
m_dot = p("m_dot", 0.212, "kg/s", "Air mass flow (HX35 map)")
rho3 = p("rho3", 2.118, "kg/m3", "Compressor delivery air density")
V_an_min = p("V_an_min", 25.0, "m/s", "Annulus air speed target, low")
V_an_max = p("V_an_max", 30.0, "m/s", "Annulus air speed target, high")
far = p("far", 0.01213, "-", "Overall fuel/air ratio")
P3 = p("P3", 243180.0, "Pa", "Combustor inlet pressure")
T3 = p("T3", 400.0, "K", "Combustor inlet temperature")
liner_dp = p("liner_dp", 0.04, "-", "Liner pressure drop / P3")
R_gas = p("R_gas", 288.0, "J/kg K", "Hot gas constant (products ~ air)")
air_prim = p("air_prim", 0.189, "-", "Cumulative air fraction in the primary zone")
air_sec = p("air_sec", 0.315, "-", "Cumulative air fraction in the secondary zone")
T_prim = p("T_prim", 2337.0, "K", "Primary zone temperature (Cantera)")
T_sec = p("T_sec", 1777.0, "K", "Secondary zone temperature (Cantera)")
T4 = p("T4", 880.0, "K", "Turbine inlet temperature target")
L_prim = p("L_prim", 70.0, "mm", "Primary zone length (x = 0 .. L_prim)")
L_sec = p("L_sec", 60.0, "mm", "Secondary zone length")
L_dil = p("L_dil", liner_L - L_prim - L_sec, "mm", "Dilution zone length (rest of the liner)", "derived")
tau_prim = p("tau_prim", 2.0e-3, "s", "Primary residence time target (~20x propane chemistry)")
tau_sec = p("tau_sec", 1.5e-3, "s", "Secondary residence time target (CO burnout)")
dil_LD_min = p("dil_LD_min", 1.5, "-", "Dilution mixing length, min (in liner diameters)")
dil_LD_max = p("dil_LD_max", 1.8, "-", "Dilution mixing length, max (in liner diameters)")
alpha_304 = p("alpha_304", 18e-6, "1/K", "304 thermal expansion")
dT_liner = p("dT_liner", 700.0, "K", "Liner temperature rise, hot vs. assembly")
T_amb = p("T_amb", 293.15, "K", "Assembly temperature")
T_sleeve_min = p("T_sleeve_min", 400.0, "K", "Slip sleeve temperature, low end (~T3, cooled by annulus air)", "assumed")
T_sleeve_max = p("T_sleeve_max", 700.0, "K", "Slip sleeve temperature, high end (heated by the liner it touches)", "assumed")
radial_hot_min = p("radial_hot_min", 0.1, "mm", "WARN if hot radial slip clearance in the T_sleeve_min case falls below this")
leak_frac_max = p("leak_frac_max", 0.05, "-", "Reference level for slip-sleeve leak area / liner hole area (reported, INFO only)")
clock_min_sep = p("clock_min_sep", 2.5, "deg", "WARN if two non-cooling hole rows are closer than this in angle")
rho_304 = p("rho_304", 7.93e-6, "kg/mm3", "304 density (mass estimate)")

p("material", "AISI 304 body; E308L filler; graphite gaskets", "", "Materials (note only)", "note")
p("bolt_note", "Match-drill the 3 bolt holes through liner + ring + spacers in one setup", "", "Shop note", "note")
p("slip_note", "Apply anti-seize / graphite paste to liner OD in the sleeve (galling)", "", "Shop note", "note")

# =============================================================================
# Geometry helpers
# =============================================================================
X_AXIS = cq.Vector(1, 0, 0)
ORIGIN = cq.Vector(0, 0, 0)


def angles(n, phase_deg):
    return [phase_deg + 360.0 * i / n for i in range(n)]


def yz(r, theta_deg):
    """Local (Y, Z) of a point at radius r, angle theta (from +Z toward +Y)."""
    t = math.radians(theta_deg)
    return (r * math.sin(t), r * math.cos(t))


def tube(x0, length, od, id_):
    """Round tube on the X axis from x0 to x0 + length."""
    return (cq.Workplane("YZ").workplane(offset=x0)
            .circle(od / 2).circle(id_ / 2).extrude(length))


def rod(x0, length, d):
    return cq.Workplane("YZ").workplane(offset=x0).circle(d / 2).extrude(length)


def to_radial(shape, x, theta_deg):
    """Move a shape built along +Z (radial direction at theta = 0) to angle
    theta about the X axis and to axial position x."""
    return shape.rotate(ORIGIN, X_AXIS, -theta_deg).translate(cq.Vector(x, 0, 0))


def radial_cyl(x, theta_deg, d, r0, r1):
    """Radial cylinder (drill) of diameter d from radius r0 to r1."""
    c = cq.Solid.makeCylinder(d / 2, r1 - r0, cq.Vector(0, 0, r0), cq.Vector(0, 0, 1))
    return to_radial(c, x, theta_deg)


def hex_prism(af, h, z0):
    """Hex prism along +Z (across flats af, height h) starting at z = z0."""
    return (cq.Workplane("XY").workplane(offset=z0)
            .polygon(6, af * 2 / math.sqrt(3)).extrude(h).val())


def cut_all(wp, cutters):
    return wp.cut(cq.Workplane().add(cq.Compound.makeCompound(cutters)))


# Hole table: one record per drilled hole in the liner. Used both to cut the
# liner and by the checks, so geometry and checks cannot drift apart.
ROWS = [
    # name,         zone,        n,         d,      x,      phase
    ("cooling",     "cooling",   n_cool,    d_cool, x_cool, phase_cool),
    ("primary",     "primary",   n_prim,    d_prim, x_prim, phase_prim),
    ("secondary",   "secondary", n_sec,     d_sec,  x_sec,  phase_sec),
    ("dilution_1",  "dilution",  n_dil_row, d_dil,  x_dil1, phase_dil1),
    ("dilution_2",  "dilution",  n_dil_row, d_dil,  x_dil2, phase_dil2),
]
HOLES = [dict(row=r, zone=z, d=d, x=x, theta=t)
         for (r, z, n, d, x, ph) in ROWS for t in angles(n, ph)]
BOLT_ANGLES = angles(n_bolt, phase_bolt)

# Bolt head sits flat inside the ring bore: lift its bearing face so the hex
# corners clear the curved bore instead of digging into it.
bolt_seat_r = math.sqrt((ring_ID / 2) ** 2 - (bolt_head_AF / math.sqrt(3)) ** 2)
_bolt_need = liner_OD / 2 - bolt_seat_r + nut_m + 1.5   # + ~2 threads past the nut
bolt_len = p("bolt_len", min(l for l in (8, 10, 12, 14, 16, 20, 25) if l >= _bolt_need),
             "mm", "M4 bolt length (next standard length that shows ~2 threads past the nut)", "derived")


# =============================================================================
# Parts (all positioned in the global frame)
# =============================================================================
def flange_disc(x0, t, bore):
    """EN 1092-1 PN6 DN125 disc from x0 to x0 + t with the 8 bolt holes.
    Holes straddle the vertical so the stub sits between two bolts."""
    disc = tube(x0, t, flange_OD, bore) if bore else rod(x0, t, flange_OD)
    pts = [yz(flange_BCD / 2, a) for a in angles(n_flange_bolt, flange_phase)]
    holes = (cq.Workplane("YZ").workplane(offset=x0 - 1)
             .pushPoints(pts).circle(flange_hole_d / 2).extrude(t + 2))
    return disc.cut(holes)


def blind_flange():
    """1. Blind flange = head plate. Inner face at x = 0; the head ring is
    welded to it, so liner + ring come out together through the DN125 pipe."""
    return flange_disc(-flange_t, flange_t, None)


def flange_gasket():
    """2. Full-face graphite gasket between blind and slip-on flange."""
    return flange_disc(0, fl_gasket_t, flange_bore)


def slipon_flange():
    """3. EN 1092-1 type 01 slip-on flange, welded to the head pipe."""
    return flange_disc(fl_gasket_t, flange_t, flange_bore)


def _stub_bore():
    return (cq.Workplane("XY").center(stub_x, 0).circle(stub_ID / 2)
            .extrude(head_OD / 2 + stub_protrude))


def head_pipe():
    """4. DN125 head pipe: slip-on flange to step plate, with the stub hole."""
    x0 = fl_gasket_t + slipon_setback
    return tube(x0, head_L - x0, head_OD, head_ID).cut(_stub_bore())


def inlet_stub():
    """5. Set-on branch: the stub end is saddle-cut to the head pipe OD."""
    stub = (cq.Workplane("XY").center(stub_x, 0)
            .circle(stub_OD / 2).circle(stub_ID / 2)
            .extrude(head_OD / 2 + stub_protrude))
    return stub.cut(rod(0, head_L + 10, head_OD))


def step_plate():
    """6. Flat annular plate DN125 -> 108 housing. A flat plate instead of a
    cone because a non-CNC shop can cut it from plate and weld it square."""
    return tube(head_L, step_t, head_OD, housing_OD)


def housing():
    """7. 108 x 2.5 housing tube: passes through the step plate, ends on the
    end plate."""
    return tube(head_L, x_end_in - head_L, housing_OD, housing_ID)


def end_plate():
    """8. Turbine-end plate with the gas outlet bore."""
    return tube(x_end_in, end_t, housing_OD, end_hole_d)


def head_ring():
    """9. 63.5 head ring, welded to the head plate. Its OD and the liner ID
    form the film slot; bolt holes match-drilled with the liner."""
    ring = tube(0, ring_L, ring_OD, ring_ID)
    return cut_all(ring, [radial_cyl(x_bolt, t, bolt_hole_d, ring_ID / 2 - 1, ring_OD / 2 + 1)
                          for t in BOLT_ANGLES])


def tube_gasket():
    """10. Graphite ring under the liner end: closes the film slot at the
    head end so all film air comes through the cooling holes."""
    return tube(0, tube_gasket_t, liner_OD, ring_OD)


def main_liner():
    """11. Main liner with all hole rows and the bolt holes."""
    liner = tube(liner_start, main_L, liner_OD, liner_ID)
    r0, r1 = liner_ID / 2 - 1, liner_OD / 2 + 1   # drill through the liner wall only
    cutters = [radial_cyl(h["x"], h["theta"], h["d"], r0, r1) for h in HOLES]
    cutters += [radial_cyl(x_bolt, t, bolt_hole_d, r0, r1) for t in BOLT_ANGLES]
    return cut_all(liner, cutters)


def spacers():
    """12. Distance pieces in the film slot, one per bolt: keep the slot
    concentric and take the bolt clamp load. Saddle-shaped to sit on both
    tubes (in practice: short tube or washer stack ground to fit)."""
    slot = tube(x_bolt - spacer_OD, 2 * spacer_OD, liner_ID, ring_OD).val()
    out = []
    for t in BOLT_ANGLES:
        blank = radial_cyl(x_bolt, t, spacer_OD, ring_OD / 2 - 1, liner_ID / 2 + 1)
        s = blank.intersect(slot).cut(
            radial_cyl(x_bolt, t, bolt_hole_d, ring_OD / 2 - 1, liner_ID / 2 + 1))
        out.append(s)
    return cq.Workplane().add(cq.Compound.makeCompound(out))


def ring_bolts():
    """13. M4 hex bolts (head inside the ring) with nuts outside the liner."""
    out = []
    for t in BOLT_ANGLES:
        head = hex_prism(bolt_head_AF, bolt_head_k, bolt_seat_r - bolt_head_k)
        shank = cq.Solid.makeCylinder(bolt_d / 2, bolt_len, cq.Vector(0, 0, bolt_seat_r))
        nut = (cq.Workplane().add(hex_prism(nut_AF, nut_m, liner_OD / 2))
               .cut(cq.Workplane().add(cq.Solid.makeCylinder(
                   bolt_d / 2, nut_m + 2, cq.Vector(0, 0, liner_OD / 2 - 1)))).val())
        out.append(to_radial(head.fuse(shank).clean(), x_bolt, t))
        out.append(to_radial(nut, x_bolt, t))
    return cq.Workplane().add(cq.Compound.makeCompound(out))


def slip_sleeve():
    """14. Sleeve on the end plate; the liner end slides in it so the liner
    is fixed at the head end only and can grow freely."""
    return tube(x_end_in - sleeve_L, sleeve_L, sleeve_OD, sleeve_ID)


def fluid_volume(parts):
    """15. Air annulus + gas inside the liner as one solid (for CFD).
    The cavity is bounded by the inner surfaces of the pressure casing
    (gasket/flange bore, head pipe, stub, housing, end plate outlet); the
    internal parts are then subtracted. Pipe-in-flange-bore clearance is
    treated as closed by the internal fillet weld."""
    x_pipe = fl_gasket_t + slipon_setback
    stub_top = head_OD / 2 + stub_protrude
    cavity = (rod(0, x_pipe, flange_bore)
              .union(rod(x_pipe, head_L - x_pipe, head_ID))
              .union(cq.Workplane("XY").center(stub_x, 0).circle(stub_ID / 2).extrude(stub_top))
              .union(rod(head_L, x_end_in - head_L, housing_ID))
              .union(rod(x_end_in, end_t, end_hole_d)))
    for name in ("09_head_ring", "10_tube_gasket", "11_main_liner",
                 "12_spacers", "13_ring_bolts", "14_slip_sleeve"):
        cavity = cavity.cut(parts[name])
    return cavity


PART_BUILDERS = [
    ("01_blind_flange", blind_flange),
    ("02_flange_gasket", flange_gasket),
    ("03_slipon_flange", slipon_flange),
    ("04_head_pipe", head_pipe),
    ("05_inlet_stub", inlet_stub),
    ("06_step_plate", step_plate),
    ("07_housing", housing),
    ("08_end_plate", end_plate),
    ("09_head_ring", head_ring),
    ("10_tube_gasket", tube_gasket),
    ("11_main_liner", main_liner),
    ("12_spacers", spacers),
    ("13_ring_bolts", ring_bolts),
    ("14_slip_sleeve", slip_sleeve),
]
GASKETS = {"02_flange_gasket", "10_tube_gasket"}


# =============================================================================
# Checks
# =============================================================================
def run_checks(parts, fluid):
    rows = []

    def chk(ok, name, detail, warn_level="WARN"):
        status = ok if isinstance(ok, str) else ("PASS" if ok else warn_level)
        rows.append((status, name, detail))

    eps = 1e-9

    # --- film slot
    chk(liner_ID - ring_OD >= 0, "Ring fits inside liner",
        f"liner_ID - ring_OD = {liner_ID - ring_OD:.2f} mm (>= 0)")
    chk(3.0 <= film_gap <= 4.5, "Film gap 3.0-4.5 mm", f"film_gap = {film_gap:.2f} mm")

    # --- annulus
    gap = (housing_ID - liner_OD) / 2
    A_an = math.pi / 4 * (housing_ID ** 2 - liner_OD ** 2)
    V_an = m_dot / (rho3 * A_an * 1e-6)
    chk(V_an_min <= V_an <= V_an_max, "Annulus air speed 25-30 m/s",
        f"gap {gap:.2f} mm, area {A_an:.1f} mm2, V = {V_an:.2f} m/s, "
        f"q = {0.5 * rho3 * V_an ** 2:.0f} Pa")
    A_head = math.pi / 4 * (head_ID ** 2 - liner_OD ** 2)
    chk("INFO", "Head-section annulus (plenum)",
        f"area {A_head:.0f} mm2, V = {m_dot / (rho3 * A_head * 1e-6):.1f} m/s")

    # --- hole areas
    targets = {"cooling": A_cool_target, "primary": A_prim_target,
               "secondary": A_sec_target, "dilution": A_dil_target}
    total = 0.0
    for zone, tgt in targets.items():
        a = sum(math.pi * h["d"] ** 2 / 4 for h in HOLES if h["zone"] == zone)
        n = sum(1 for h in HOLES if h["zone"] == zone)
        total += a
        dev = a / tgt - 1
        chk(abs(dev) <= area_tol, f"Hole area {zone}",
            f"{n} holes = {a:.1f} mm2 vs target {tgt:.1f} ({dev * 100:+.2f} %, tol {area_tol * 100:.0f} %)")
    tgt_total = sum(targets.values())
    chk(abs(total / tgt_total - 1) <= area_tol, "Hole area total",
        f"{total:.1f} mm2 vs {tgt_total:.1f} mm2 ({(total / tgt_total - 1) * 100:+.2f} %)")

    # --- hole overlaps (unrolled liner OD surface, the worst case for ligaments)
    r_s = liner_OD / 2
    feats = [dict(h, kind="hole") for h in HOLES]
    feats += [dict(row="bolt", zone="bolt", d=bolt_hole_d, x=x_bolt, theta=t, kind="bolt") for t in BOLT_ANGLES]
    # spacer footprint on the liner ID: no hole may open onto it
    feats += [dict(row="spacer", zone="spacer", d=spacer_OD, x=x_bolt, theta=t, kind="spacer") for t in BOLT_ANGLES]
    worst = (float("inf"), None, None)
    for a, b in combinations(feats, 2):
        if a["kind"] != "hole" and b["kind"] != "hole":
            continue   # bolt vs. its own spacer is intended
        dth = math.radians((a["theta"] - b["theta"] + 180) % 360 - 180)
        dist = math.hypot(a["x"] - b["x"], r_s * dth)
        lig = dist - (a["d"] + b["d"]) / 2
        if lig < worst[0]:
            worst = (lig, a, b)
    lig, a, b = worst
    chk(lig > 0, "No hole overlaps hole / bolt / spacer",
        f"min ligament {lig:.2f} mm ({a['row']} @ {a['theta']:.1f} deg vs {b['row']} @ {b['theta']:.1f} deg)")

    # --- clocking: minimum angular separation between every pair of rows.
    # A jet directly behind an upstream one sits in its wake. The cooling row is
    # exempt (INFO): it feeds the film slot under the ring and leaves no wake.
    for (ra, za, na, _, _, pa), (rb, zb, nb, _, _, pb) in combinations(ROWS, 2):
        sep = min(abs((ta - tb + 180) % 360 - 180)
                  for ta in angles(na, pa) for tb in angles(nb, pb))
        if "cooling" in (za, zb):
            chk("INFO", f"Clocking {ra} / {rb}", f"min angular separation {sep:.3f} deg (cooling row exempt)")
        else:
            chk(sep >= clock_min_sep, f"Clocking {ra} / {rb}",
                f"min angular separation {sep:.3f} deg (WARN < {clock_min_sep})")

    # --- ring overlap: only the cooling row may sit over the ring
    zone_holes = [r for r in ROWS if r[1] != "cooling"]
    first_edge = min(x - d / 2 for (_, _, _, d, x, _) in zone_holes)
    chk(first_edge > ring_L, "Zone holes clear of ring overlap",
        f"first zone-hole edge x = {first_edge:.2f} > ring end {ring_L:.1f} "
        f"(margin {first_edge - ring_L:.2f} mm)")
    c_lo, c_hi = x_cool - d_cool / 2, x_cool + d_cool / 2
    chk(liner_start < c_lo and c_hi < ring_L, "Cooling row inside film path",
        f"cooling holes x {c_lo:.1f}-{c_hi:.1f} within overlap {liner_start:.1f}-{ring_L:.1f}; "
        f"drilled through liner only (ring not cut); film runs {ring_L - c_hi:.1f} mm to the slot exit")
    b_hi = x_bolt + max(bolt_hole_d, spacer_OD) / 2
    chk(liner_start < x_bolt - spacer_OD / 2 and b_hi < c_lo, "Bolt row inside overlap, upstream of cooling",
        f"bolt/spacer x {x_bolt - spacer_OD / 2:.1f}-{b_hi:.1f}, cooling starts {c_lo:.1f}")

    # --- every row inside the liner and clear of ring end and slip sleeve
    sleeve_start = x_end_in - sleeve_L
    for (name, zone, n, d, x, ph) in ROWS:
        lo, hi = x - d / 2, x + d / 2
        if zone == "cooling":
            ok = liner_start < lo and hi < ring_L
            where = f"in overlap {liner_start:.0f}-{ring_L:.0f}"
        else:
            ok = ring_L < lo and hi < sleeve_start
            where = f"between ring end {ring_L:.0f} and sleeve start {sleeve_start:.0f}"
        chk(ok, f"Row x-position {name}", f"x {lo:.2f}-{hi:.2f} {where}")

    # --- rows sit in their zones
    bounds = {"primary": (0, L_prim), "secondary": (L_prim, L_prim + L_sec),
              "dilution": (L_prim + L_sec, liner_L)}
    for (name, zone, n, d, x, ph) in ROWS:
        if zone in bounds:
            lo, hi = bounds[zone]
            chk(lo <= x - d / 2 and x + d / 2 <= hi, f"Row {name} inside {zone} zone",
                f"x = {x:.1f} within {lo:.0f}-{hi:.0f}")

    # --- head section layout
    flange_back = fl_gasket_t + flange_t
    stub_edge = stub_x - stub_OD / 2
    chk(stub_edge - flange_back >= hw_clear - eps, "Stub clear of flange hardware",
        f"stub edge x = {stub_edge:.2f}, flange back x = {flange_back:.2f}, clearance "
        f"{stub_edge - flange_back:.2f} >= hw_clear {hw_clear:.1f}")
    chk(abs(stub_x - (fl_gasket_t + flange_t + hw_clear + stub_OD / 2)) < eps, "stub_x formula",
        f"stub_x = {fl_gasket_t} + {flange_t} + {hw_clear} + {stub_OD / 2} = {stub_x:.2f}")
    chk(abs(head_L - (stub_x + stub_OD / 2 + 10)) < eps, "head_L formula",
        f"head_L = {stub_x:.2f} + {stub_OD / 2} + 10 = {head_L:.2f}")
    chk(head_L - (stub_x + stub_OD / 2) >= 10 - eps, "Stub clear of step plate",
        f"{head_L - (stub_x + stub_OD / 2):.2f} mm weld room")

    # --- residence time: L = tau * V, then tau at the chosen lengths
    A_l = math.pi / 4 * (liner_ID * 1e-3) ** 2
    P_l = P3 * (1 - liner_dp)
    m_f = far * m_dot

    def speed(m, T):
        return m / (P_l / (R_gas * T) * A_l)

    V_p = speed(air_prim * m_dot + m_f, T_prim)
    V_s = speed(air_sec * m_dot + m_f, T_sec)
    V_d_out = speed(m_dot + m_f, T4)
    V_d = (V_s + V_d_out) / 2       # dilution gas accelerates from V_s to exit speed
    for zone, V, tau, L in (("primary", V_p, tau_prim, L_prim), ("secondary", V_s, tau_sec, L_sec)):
        L_req = tau * V * 1e3
        chk(L >= L_req, f"Residence length {zone}",
            f"V = {V:.1f} m/s, L_req = tau*V = {tau * 1e3:.1f} ms x V = {L_req:.1f} mm <= {L:.0f} mm "
            f"(tau = {L / V:.2f} ms)")
    LD = L_dil / liner_ID
    chk(dil_LD_min <= LD <= dil_LD_max, "Dilution mixing length",
        f"L_dil = {L_dil:.0f} mm = {LD:.2f} D (target {dil_LD_min}-{dil_LD_max} D), "
        f"V {V_s:.1f} -> {V_d_out:.1f} m/s, tau = {L_dil / V_d:.2f} ms")
    tau_tot = L_prim / V_p + L_sec / V_s + L_dil / V_d
    chk("INFO", "Residence time at liner_L",
        f"tau = {tau_tot:.2f} ms over liner_L = {liner_L:.0f} mm (plug flow, film neglected)")

    # --- slip joint thermal growth
    # Liner growth against sleeve growth for both ends of the assumed sleeve
    # temperature range. A cool sleeve (annulus air) is the tight case.
    dD_liner = liner_OD * alpha_304 * dT_liner
    A_holes = sum(math.pi * h["d"] ** 2 / 4 for h in HOLES)
    A_leak = math.pi * (sleeve_ID ** 2 - liner_OD ** 2) / 4
    chk("INFO", "Slip sleeve leak area, cold",
        f"{A_leak:.1f} mm2 = {A_leak / A_holes * 100:.2f} % of liner hole area {A_holes:.1f} mm2 "
        f"(reference {leak_frac_max * 100:.0f} %)")
    for T_s in (T_sleeve_min, T_sleeve_max):
        dD_sleeve = sleeve_ID * alpha_304 * (T_s - T_amb)
        rad_hot = slip_radial - (dD_liner - dD_sleeve) / 2
        detail = (f"liner growth {dD_liner:.3f} mm, sleeve growth {dD_sleeve:.3f} mm (dia.); "
                  f"hot radial clearance = {slip_radial} - ({dD_liner:.3f} - {dD_sleeve:.3f})/2 = {rad_hot:.3f} mm")
        if T_s == T_sleeve_min:
            chk(rad_hot >= radial_hot_min, f"Slip hot radial clearance, sleeve {T_s:.0f} K",
                detail + f" (WARN < {radial_hot_min})")
        else:
            chk("INFO", f"Slip hot radial clearance, sleeve {T_s:.0f} K", detail)
        D_s, D_l = sleeve_ID + dD_sleeve, liner_OD + dD_liner
        A_hot = math.pi * (D_s ** 2 - D_l ** 2) / 4
        chk("INFO", f"Slip sleeve leak area, hot, sleeve {T_s:.0f} K",
            f"{A_hot:.1f} mm2 = {A_hot / A_holes * 100:.2f} % of liner hole area "
            f"(reference {leak_frac_max * 100:.0f} %)")
    dL_liner = main_L * alpha_304 * dT_liner
    ax_hot = slip_axial - dL_liner
    chk(ax_hot > 0, "Slip joint hot axial clearance",
        f"liner hot length growth {dL_liner:.2f} mm; hot axial clearance = {slip_axial} - "
        f"{dL_liner:.2f} = {ax_hot:.2f} mm (casing growth not credited); cold engagement "
        f"{sleeve_L - slip_axial:.1f} mm")

    # --- inlet stub and head-section spacing
    A_stub = math.pi * (stub_ID * 1e-3) ** 2 / 4
    V_stub = m_dot / (rho3 * A_stub)
    chk("INFO", "Inlet stub velocity",
        f"V = m_dot / (rho3 * pi * stub_ID^2 / 4) = {V_stub:.1f} m/s (stub_ID {stub_ID:.1f} mm), "
        f"q = {0.5 * rho3 * V_stub ** 2:.0f} Pa")
    # The stub jet's dynamic head is assumed fully lost when it dumps into the plenum.
    q_dump = 0.5 * rho3 * V_stub ** 2
    chk("INFO", "Inlet dump loss", f"0.5 * rho3 * V_stub^2 = {q_dump:.0f} Pa = {q_dump / P3 * 100:.2f} % of P3")
    chk("INFO", "Total combustor pressure loss",
        f"liner_dp {liner_dp * 100:.1f} % + dump {q_dump / P3 * 100:.2f} % = {(liner_dp + q_dump / P3) * 100:.2f} % of P3")
    chk("INFO", "Step plate to secondary row",
        f"x_sec - (head_L + step_t) = {x_sec:.1f} - {head_L + step_t:.1f} = {x_sec - head_L - step_t:.1f} mm")
    chk("INFO", "Stub far edge to step plate",
        f"head_L - (stub_x + stub_OD/2) = {head_L:.1f} - {stub_x + stub_OD / 2:.1f} = "
        f"{head_L - stub_x - stub_OD / 2:.1f} mm")

    # --- solid model sanity
    worst_int = (0.0, "")
    names = list(parts)
    for a, b in combinations(names, 2):
        sa, sb = parts[a].val(), parts[b].val()
        ba, bb = sa.BoundingBox(), sb.BoundingBox()
        if (ba.xmax < bb.xmin or bb.xmax < ba.xmin or ba.ymax < bb.ymin or
                bb.ymax < ba.ymin or ba.zmax < bb.zmin or bb.zmax < ba.zmin):
            continue   # boxes apart: skip the expensive boolean
        v = sa.intersect(sb).Volume()
        if v > worst_int[0]:
            worst_int = (v, f"{a} / {b}")
    chk(worst_int[0] < 0.01, "No part interference",
        f"max overlap {worst_int[0]:.4f} mm3 {worst_int[1]}")
    nsol = len(fluid.solids().vals())
    chk(nsol == 1, "Fluid volume is one solid",
        f"{nsol} solid(s), {fluid.val().Volume() / 1e3:.1f} cm3")
    return rows


# =============================================================================
# Outputs
# =============================================================================
COLORS = {
    "01_blind_flange": "#8a96a3", "02_flange_gasket": "#2b2b2b", "03_slipon_flange": "#8a96a3",
    "04_head_pipe": "#a9b4bf", "05_inlet_stub": "#c3ccd5", "06_step_plate": "#7d8995",
    "07_housing": "#a9b4bf", "08_end_plate": "#7d8995", "09_head_ring": "#d08a3c",
    "10_tube_gasket": "#2b2b2b", "11_main_liner": "#c0503a", "12_spacers": "#3c7dd0",
    "13_ring_bolts": "#3aa05a", "14_slip_sleeve": "#8f5cc0",
}


def write_half_section(parts, fluid, path_svg):
    """2D section of the assembly in the plane through the engine axis and the
    stub axis (XZ, y = 0), viewed from +Y. Faces filled per part, fluid in
    light blue. Drawn directly from the solids so it always matches the STEP."""
    plane = cq.Face.makePlane(4000, 4000, (0, 0, 0), (0, 1, 0))
    layers = [("15_fluid_volume", fluid, "#dcecf8")] + \
             [(n, w, COLORS[n]) for n, w in parts.items()]
    xmin, xmax = -flange_t - 10, x_end_in + end_t + 10
    zmax = max(head_OD / 2 + stub_protrude, flange_OD / 2) + 25
    zmin = -flange_OD / 2 - 10
    s = 3.0                                   # px per mm
    W, H = (xmax - xmin) * s, (zmax - zmin) * s
    legend_w = 260
    X = lambda x: (x - xmin) * s
    Z = lambda z: (zmax - z) * s
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W + legend_w:.0f}" height="{H + 75:.0f}" '
           f'viewBox="0 0 {W + legend_w:.0f} {H + 75:.0f}" font-family="Helvetica,Arial,sans-serif">',
           f'<rect width="100%" height="100%" fill="white"/>']
    for name, wp, color in layers:
        sec = wp.val().intersect(plane)
        g = [f'<g id="{name}">']
        for f in sec.Faces():
            verts, tris = f.tessellate(0.05, 0.1)
            for a, b, c in tris:
                pts = " ".join(f"{X(verts[i].x):.2f},{Z(verts[i].z):.2f}" for i in (a, b, c))
                g.append(f'<polygon points="{pts}" fill="{color}" stroke="{color}" stroke-width="0.8"/>')
        for e in sec.Edges():
            n = 2 if e.geomType() == "LINE" else 24
            pts = e.positions([i / (n - 1) for i in range(n)])
            d = " ".join(f"{X(v.x):.2f},{Z(v.z):.2f}" for v in pts)
            stroke = "#6a9cc4" if name.startswith("15") else "#1a1a1a"
            g.append(f'<polyline points="{d}" fill="none" stroke="{stroke}" stroke-width="0.6"/>')
        g.append("</g>")
        out += g
    # axis and station marks
    out.append(f'<line x1="{X(xmin):.1f}" y1="{Z(0):.1f}" x2="{X(xmax):.1f}" y2="{Z(0):.1f}" '
               f'stroke="#444" stroke-width="0.8" stroke-dasharray="12,3,2,3"/>')
    stations = [(0, "0"), (ring_L, "ring end"), (x_prim, "prim"), (x_sec, "sec"),
                (x_dil1, "dil 1"), (x_dil2, "dil 2"), (head_L, "step"), (liner_L, "liner end")]
    yb = Z(zmin) + 4
    for i, (x, lab) in enumerate(stations):
        ty = yb + 18 + 13 * (i % 2)   # alternate two label rows so close stations stay readable
        out.append(f'<line x1="{X(x):.1f}" y1="{Z(-liner_OD / 2 - 4):.1f}" x2="{X(x):.1f}" y2="{ty - 10:.1f}" '
                   f'stroke="#999" stroke-width="0.5" stroke-dasharray="3,2"/>')
        out.append(f'<text x="{X(x):.1f}" y="{ty:.1f}" font-size="10" text-anchor="middle">'
                   f'{lab} x={x:g}</text>')
    out.append(f'<text x="{X(stub_x):.1f}" y="{Z(head_OD / 2 + stub_protrude) - 4:.1f}" font-size="11" '
               f'text-anchor="middle">inlet stub x={stub_x:g}</text>')
    out.append(f'<text x="10" y="{H + 68:.0f}" font-size="12">Half-section through engine axis (X) and inlet '
               f'stub axis (Z), plane y = 0, viewed from +Y. Dimensions in mm.</text>')
    # legend
    lx, ly = W + 20, 30
    out.append(f'<text x="{lx}" y="{ly - 10}" font-size="13" font-weight="bold">Parts</text>')
    for i, (name, _, color) in enumerate(layers[1:] + layers[:1]):
        y = ly + i * 20
        out.append(f'<rect x="{lx}" y="{y}" width="14" height="14" fill="{color}" stroke="#1a1a1a" stroke-width="0.5"/>')
        out.append(f'<text x="{lx + 22}" y="{y + 11}" font-size="11">{name}</text>')
    out.append("</svg>")
    path_svg.write_text("\n".join(out))


def write_parameters(path):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "value", "unit", "kind", "description"])
        for name, (value, unit, desc, kind) in PARAMS.items():
            v = f"{value:.6g}" if isinstance(value, float) else value
            w.writerow([name, v, unit, kind, desc])


def format_checks(rows, parts, fluid):
    w1 = max(len(r[1]) for r in rows)
    lines = [f"{'STATUS':6}  {'CHECK':{w1}}  DETAIL", "-" * (w1 + 90)]
    lines += [f"{s:6}  {n:{w1}}  {d}" for s, n, d in rows]
    n_warn = sum(1 for r in rows if r[0] == "WARN")
    lines += ["", f"{sum(1 for r in rows if r[0] == 'PASS')} PASS, {n_warn} WARN, "
                  f"{sum(1 for r in rows if r[0] == 'INFO')} INFO", "",
              "Part summary (mass: 304 density; gaskets are graphite, mass not meaningful)",
              f"{'part':20} {'volume cm3':>11} {'mass kg':>8}"]
    for name, wp in parts.items():
        v = wp.val().Volume()
        m = "-" if name in GASKETS else f"{v * rho_304:.3f}"
        lines.append(f"{name:20} {v / 1e3:11.2f} {m:>8}")
    lines.append(f"{'fluid_volume':20} {fluid.val().Volume() / 1e3:11.2f} {'-':>8}")
    return "\n".join(lines)


def main():
    OUT.mkdir(exist_ok=True)
    parts = {}
    for name, build in PART_BUILDERS:
        parts[name] = build()
        print(f"built {name}")
    fluid = fluid_volume(parts)
    print("built fluid_volume")

    for name, wp in parts.items():
        cq.exporters.export(wp, str(OUT / f"{name}.STEP"))
    cq.exporters.export(fluid, str(OUT / "fluid_volume.STEP"))

    assy = cq.Assembly(name="combustor")
    for name, wp in parts.items():
        assy.add(wp, name=name, color=cq.Color(COLORS[name]))
    assy.export(str(OUT / "assembly.STEP"))

    write_parameters(OUT / "parameters.csv")
    svg = OUT / "half_section.svg"
    write_half_section(parts, fluid, svg)
    try:
        import cairosvg
        cairosvg.svg2png(url=str(svg), write_to=str(OUT / "half_section.png"))
    except ImportError:
        print("cairosvg not installed: half_section.png skipped")

    from verify_step import verify_outputs   # reads the exported STEP files back
    verify = verify_outputs(OUT, ROWS, BOLT_ANGLES, bolt_hole_d, x_bolt, parts)
    (OUT / "verify.txt").write_text(verify + "\n")

    report = format_checks(run_checks(parts, fluid), parts, fluid)
    (OUT / "checks.txt").write_text(report + "\n")
    print()
    print(report)
    print()
    print(verify)
    print("\nGenerated files:")
    for f in sorted(OUT.iterdir()):
        print(f"  output/{f.name:28} {f.stat().st_size / 1024:8.1f} kB")


if __name__ == "__main__":
    main()
