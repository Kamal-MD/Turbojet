Project description: homemade turbojet made from a salvaged turbocharger, propane fueled, ground test only\
Current Status: Found Turbocharger, broken compressor seal ring but in good condition. Combustion chamber geometry complete and verified (scripted CAD model in `Combustor_Design/`); next step is cold-flow CFD\
Current Reference: HX35 Turbo\
Key points:

## Design point summary (at full power)

### Inputs

| Parameter | Value | Source |
|---|---|---|
| Fuel | Propane (C₃H₈) | Design decision |
| Compressor pressure ratio | 2.4 | HX35 compressor map |
| Air mass flow | 0.212 kg/s (28 lb/min) | HX35 compressor map |
| Compressor inlet, T2 / P2 *(new)* | 293.15 K / 101,325 Pa | Ambient assumption |
| Compressor efficiency *(new)* | 0.78 | Back-calculated to give T3 ≈ 400 K, not measured |

### Stations

| Station | Quantity | Value |
|---|---|---|
| 3 (combustor inlet) | T3 | 400 K (≈ 127 °C) |
| 3 (combustor inlet) | P3 *(new)* | 243,180 Pa (2.43 bar) |
| 4 (turbine inlet) | T4 target | 880 K (606.7 °C) |

### Combustion

| Quantity | Value | Note |
|---|---|---|
| FAR stoichiometric | 0.0642 | Property of propane |
| φ_overall | 0.189 | Exact value; 0.19 when rounded |
| f (overall FAR) *(new)* | 0.01213 kg fuel / kg air | Exact, keeps the (1 + f) term |
| Propane flow *(new)* | 2.57 g/s | f × ṁ_air = 0.01213 × 0.212 kg/s |
| Combustion efficiency, η_b | 0.95 (assumed) | Unmeasured until testing |
| Average gas cp, T3 → T4 *(new)* | 1100 J/kg·K | Overall energy balance |
| Propane LHV *(new)* | 46.35 MJ/kg | |

### Zones

| Zone | φ | Cum. air | Air added (zone holes) | Cooling air along wall *(updated 2026-10-29)* | T (Cantera) | Role |
|---|---|---|---|---|---|---|
| Primary | 1.0 | 18.9% | 18.9% | 5.0% | 2337 K | Stable flame in recirculation bubble |
| Secondary | 0.6 (provisional) | 31.5% | 12.6% | 0% | 1777 K | Burn out CO |
| Dilution | 0.189 | 100% | 63.5% | 0% | 880 K target | Cool to T4, even out profile |
| **Total** | | | **95%** | **5%** | | |

Wall-cooling air *(updated 2026-10-29)*: 5% of total air, as one film along the primary-zone wall.
The film comes from a stepped overlap joint at the head end: the main liner slides over a shorter head ring,
and 28 × 2.0 mm holes feed the 3.8 mm slot between them. The film leaves the slot at x = 30 mm.
Wall-temperature estimates show the secondary and dilution walls survive without film, so the earlier
secondary (3%) and dilution (2%) cooling was dropped and returned to the dilution holes (58.5% → 63.5%).
Cooling air is assumed to stay along the wall and not take part in any zone's φ, so φ_overall, T4, and the
primary and secondary results are unchanged.

### Limits and margins *(new)*

| Quantity | Value | Note |
|---|---|---|
| Turbine wheel limit | 900 °C (1173 K) | Inconel 713C, conservative end of 900–982 °C |
| T4 margin, hand (η_b = 0.95) | 25.0% | 879.8 K |
| T4 margin, Cantera (ideal) | 22.7% | 907 K, upper bound |
| Liner material limit | ~870 °C (~1140 K) | 304 stainless |
| Lean flammability limit (propane) | φ ≈ 0.5 | Real limit at idle likely higher |

### Validation *(new)*

| Check | Result |
|---|---|
| Benchmark: propane-air, φ = 1, 298 K, 1 atm | 2265.7 K vs. 2260 ± 30 K: PASS |
| Primary zone dissociation | Hand overshoots Cantera by 214 K (So we will use Cantera value of T(pz)=2337K)|
| Equilibrium CO, primary → secondary | 1.27% → 26.7 ppm |

## Combustion chamber geometry *(new, 2026-10-29)*

Built by `Combustor_Design/combustor_model.py` (CadQuery), which exports every part as STEP, the
assembly, the air volume for CFD, and `parameters.csv` with every dimension. The script and its outputs
are the source of truth for the geometry; Fusion 360 is used to view, check and draw from the STEP files.
The script was written by Claude Code under my direction: I set the design, made every decision, and
reviewed each revision (see `log/2026-10-29.md`). Details: [`Combustor_Design/README.md`](Combustor_Design/README.md).

Axis is X, with x = 0 at the inner face of the head plate. All material 304 stainless, E308L filler, graphite gaskets.

### Construction

| Part | Size | Note |
|---|---|---|
| Main liner | 76.1 × 2.5 mm tube, 249 mm long | Butts on a graphite gasket at the head plate; fixed only at the head end |
| Head ring | 63.5 × 2.5 mm tube, 30 mm long | Welded to the head plate; forms the 3.8 mm film slot inside the liner |
| Liner retention | 3 × M4 bolts at x = 9 mm with spacers | Bolt head inside, nut outside; holes match-drilled through liner and ring |
| Housing | 108 × 2.5 mm tube | Annulus gap 13.45 mm, ~26 m/s |
| Head section | DN125 pipe (139.7 mm OD), flat 4 mm step plate down to the housing | Plenum, ~10 m/s; no cone or forming |
| Head joint | EN 1092-1 PN6 DN125 blind + slip-on flange, 8 × M16 | Bolted, so the liner comes out with the head plate |
| Inlet stub | 50 OD / 45 ID, at x = 61 mm | Fits the existing 50.5 mm ID silicone hose; ~63 m/s |
| Turbine end | Slip sleeve, 0.5 mm radial / 5 mm axial clearance | Liner grows ~3.1 mm long and ~1.0 mm in diameter when hot |

### Liner holes

| Group | Holes | x | First angle | Area vs. target | Pitch (hole diameters, on ID) |
|---|---|---|---|---|---|
| Cooling (film slot) | 28 × 2.0 mm | 18 mm | 0° | +1.0% | 4.0 |
| Primary | 10 × 6.5 mm | 62 mm | 0° | +0.9% | 3.4 |
| Secondary | 14 × 4.5 mm | 115 mm | 12.857° | +1.5% | 3.5 |
| Dilution (2 staggered rows of 7) | 14 × 10 mm | 145 / 170 mm | 2.571° / 28.286° | −0.5% | 3.2 |

Total hole area 1742 mm² (target 1741 mm²). Rows are clocked so no gas-side jet sits directly behind
another; 2.571° is the best possible separation between 10-hole and 14-hole rows.

### Checks (from `Combustor_Design/output/checks.txt`)

| Check | Result |
|---|---|
| STEP read-back of all hole rows | All match the parameters |
| Part interference (21 solids) | None; all contacts intended |
| Total combustor pressure loss | 5.72% of P3 (4% liner + 1.72% inlet dump) |
| Slip sleeve hot radial clearance | 0.095 mm (sleeve 400 K) to 0.303 mm (sleeve 700 K): WARN, set the sleeve bore from the measured liner OD |
| Sleeve leak | 6.9% of hole area cold, 1.3–4.2% hot |

### Open items

- Cold-flow CFD in SimScale on `fluid_volume.STEP`: hole feed around each row, stub jet at the primary row, step-plate separation at the secondary row
- Rerun the cycle with 5.72% combustor pressure loss (was 4%)
- Measure real tubes, pipe and flanges and update the placeholders (`flange_t`, gasket thicknesses, `head_wall`)
- Not modelled yet: exit transition to the T3 turbine flange, fuel injector, igniter
- Shop drawings, one per part, after CFD

Note: Not for flight
