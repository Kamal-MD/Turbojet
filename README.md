Project description: homemade turbojet made from a salvaged turbocharger, propane fueled, ground test only\
Current Status: Found Turbocharger, broken compressor seal ring but in good condition\
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

| Zone | φ | Cum. air | Air added | T (Cantera) | Role |
|---|---|---|---|---|---|
| Primary | 1.0 | 18.9% | 18.9% | 2337 K | Stable flame in recirculation bubble |
| Secondary | 0.6 (provisional) | 31.5% | 12.6% | 1777 K | Burn out CO |
| Dilution | 0.189 | 100% | 68.5% | 880 K target | Cool to T4, even out profile |

Wall-cooling air: not yet allocated (0%).

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
Note: Not for flight

