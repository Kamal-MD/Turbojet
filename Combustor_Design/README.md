# Combustor design: parametric CadQuery model

`combustor_model.py` builds a reference model of the propane combustor (ground test, 304 stainless,
made from standard tube, pipe and flanges) and exports STEP files to import into Fusion 360 as a
reference. Dimensions are in mm.

Frame: X is the engine axis and flow goes toward +X. x = 0 is the inner face of the head plate
(the blind flange). The inlet stub points along +Z. Hole angles are measured about X, from +Z toward +Y.

## Run

CadQuery needs Python 3.10–3.12:

```bash
cd Combustor_Design
uv venv --python 3.11 .venv          # or: python3.11 -m venv .venv
uv pip install -r requirements.txt   # or: .venv/bin/pip install -r requirements.txt
.venv/bin/python combustor_model.py
```

It runs in about 10 s, prints the check table and the generated files, and writes everything to `output/`.
Every input is set in the `PARAMETERS` block at the top of the script.

## Outputs (`output/`)

| File | What it is |
|---|---|
| `01_blind_flange.STEP` | Blind flange (the head plate), PN6 DN125, inner face at x = 0 |
| `02_flange_gasket.STEP` | Full-face graphite flange gasket |
| `03_slipon_flange.STEP` | Slip-on flange (EN 1092-1 type 01) on the head pipe |
| `04_head_pipe.STEP` | DN125 head pipe with the stub hole |
| `05_inlet_stub.STEP` | 55/50 inlet stub, saddle-cut to the head pipe (set-on branch) |
| `06_step_plate.STEP` | Flat annular plate from DN125 down to the 108 housing |
| `07_housing.STEP` | 108 × 2.5 housing tube |
| `08_end_plate.STEP` | End plate at the turbine end, with the gas outlet bore |
| `09_head_ring.STEP` | 63.5 head ring, welded to the head plate |
| `10_tube_gasket.STEP` | Graphite ring under the upstream end of the liner |
| `11_main_liner.STEP` | Main liner with all hole rows and the bolt holes |
| `12_spacers.STEP` | Three film-gap spacers, one at each bolt |
| `13_ring_bolts.STEP` | 3 × M4 hex bolts (head inside the ring) and nuts |
| `14_slip_sleeve.STEP` | Turbine-end slip sleeve on the end plate |
| `fluid_volume.STEP` | Part 15: air annulus and the gas inside the liner as one solid, for CFD |
| `assembly.STEP` | Parts 1–14 in position, named and coloured (the fluid volume is not included) |
| `half_section.svg` / `.png` | Section through the engine axis and the stub axis (plane y = 0) |
| `parameters.csv` | Every parameter with its value, unit, kind and description |
| `checks.txt` | PASS/WARN check table, plus volume and mass of each part |

The PNG is not committed: it is regenerated from the SVG whenever `cairosvg` is installed.

## Parameter kinds (`parameters.csv`)

- **input**: from the design spec. Change one only after reviewing the design.
- **derived**: calculated from the inputs.
- **assumed**: a placeholder the model needed but the spec does not give, such as the spacer OD,
  sleeve length, stub height, hole clocking or end plate thickness. Review these.
- **standard**: a value taken from a standard (fastener sizes, the flange bore). Check it against the
  supplier's catalogue.

## Assembly logic

- The head ring is welded to the blind flange. The liner slides over the ring onto the tube gasket and
  is held by 3 × M4 bolts at x = 9, with spacers in the 3.8 mm film slot. Removing the blind flange
  brings the liner out through the DN125 pipe.
- Cooling air enters through 28 × 2.0 holes at x = 18, drilled through the liner only. It runs along
  the slot and leaves as a wall film at the ring end (x = 30).
- At the turbine end the liner slides in a sleeve, with 0.5 mm radial and 5 mm axial clearance. The
  liner is therefore fixed only at the head end and can grow freely when hot.
- Not modelled: fuel injector, igniter, M16 flange bolts, welds, and the exit transition to the HX35 T3
  turbine flange.

Materials: 304 body, E308L filler, graphite gaskets.
