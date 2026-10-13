"""
chamber_calc.py - Combustion chamber design calculator for a homemade turbojet.

Station numbering:
    2 = compressor inlet
    3 = compressor exit / combustor inlet
    4 = combustor exit / turbine inlet

Edit the INPUTS block below to update the design. Everything else is computed.
Run with:  python calc/chamber_calc.py
"""

import sys
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # save plots to file, never open a window
import matplotlib.pyplot as plt
import numpy as np
import cantera as ct
from scipy.optimize import brentq

sys.stdout.reconfigure(encoding="utf-8")  # so φ and ° print on Windows consoles


# =============================================================================
# INPUTS  (all SI units; change one line to update the design)
# =============================================================================
m_air = 0.212               # kg/s, air mass flow (HX35 map: 28 lb/min)
T2 = 293.15                 # K, compressor inlet temperature (ambient)
P2 = 101325                 # Pa, compressor inlet pressure
PR = 2.4                    # -, compressor pressure ratio (HX35 map)
eta_c = 0.78                # -, compressor isentropic efficiency (gives T3 ≈ 400 K)
T4_target = 879.85          # K, turbine inlet target (606.7 °C)
T_turbine_max = 1173.15     # K, turbine material limit (900 °C)
T4_margin_min = 0.25        # -, required margin below T_turbine_max (25 %)
phi_overall_design = None   # -, rounded design value; None = use the HAND-calculated value
phi_primary = 1.0           # -, primary zone equivalence ratio
phi_secondary = 0.6         # -, secondary zone equivalence ratio (provisional)
cooling_fraction = 0.10     # -, fraction of TOTAL air used for liner wall cooling (was 0.0)
# Assumption: cooling air runs along the inner liner wall as a film and takes NO
# part in any zone's φ; it mixes into the gas only before the turbine. So the
# primary/secondary φ, cumulative air and temperatures are unchanged; only the
# dilution-hole air shrinks. φ_overall and T4 are unchanged (same total air and fuel).
LHV = 46.35e6               # J/kg, propane lower heating value
eta_b_exit = 0.95           # -, combustion efficiency for the overall energy balance
eta_b_hot = 1.0             # -, ideal efficiency for hot-zone hand calcs (comparable to Cantera)
cp_overall = 1100           # J/kg·K, average gas cp from T3 to T4 (overall balance)
cp_primary = 1300           # J/kg·K, average hot-gas cp, primary zone
cp_secondary = 1250         # J/kg·K, average hot-gas cp, secondary zone
cp_dilution = 1100          # J/kg·K, average hot-gas cp, dilution zone
MECH = "gri30.yaml"         # Cantera reaction mechanism (GRI-Mech 3.0, includes propane)
FUEL = "C3H8"               # fuel species (propane)
OXIDIZER = "O2:1, N2:3.76"  # air as a molar O2/N2 mixture
# =============================================================================


# =============================================================================
# LINER PRESSURE DROP AND HOLE SIZING INPUTS  (added; existing inputs above unchanged)
# =============================================================================
dP_liner_fraction = 0.04        # -, liner pressure drop as a fraction of P3 (design choice)
dP_sweep = [0.03, 0.04, 0.05]   # -, dP fractions to compare in the report
Cd = 0.6                        # -, discharge coefficient, sharp-edged hole value (assumed; real range ~0.5-0.8)
R_air = 287.05                  # J/kg·K, gas constant for air
d_min_mm = 2.5                  # mm, smallest practical hole (not below plate thickness)
d_max_mm = 15.0                 # mm, largest practical hole
d_step_mm = 0.5                 # mm, diameter increment available for drilling
liner_diameter_mm = None        # mm, not chosen yet; when None, skip the per-row fit check
spacing_factor = 3.0            # -, min center-to-center spacing, in hole diameters
cooling_distribution = {"primary": 0.5, "secondary": 0.3, "dilution": 0.2}
                                # -, PROVISIONAL share of cooling air along each zone's wall;
                                #    weighted toward the hottest walls; must sum to 1
Cd_cooling = 0.6                # -, discharge coefficient for cooling holes (placeholder; small holes in 2.5 mm plate may reach ~0.75)
Cd_cooling_sensitivity = 0.75   # -, second Cd for cooling holes, report only
d_cool_min_mm = 2.0             # mm, smallest cooling hole (allowed smaller than zone holes)
d_cool_max_mm = 4.0             # mm, largest cooling hole
# =============================================================================


OUT_DIR = Path(__file__).parent / "output"
gas = ct.Solution(MECH)  # one Cantera gas object, reused everywhere


def c_of(T):
    """Convert kelvin to degrees Celsius."""
    return T - 273.15


def stoich_far():
    """Stoichiometric fuel-air mass ratio for propane.

    C3H8 + 5(O2 + 3.76 N2) -> 3 CO2 + 4 H2O + 18.8 N2: one mole of fuel needs
    exactly 5 moles of O2, so FAR_st = (mass of 1 mol fuel) / (mass of 5 mol air).
    """
    M_C, M_H, M_O, M_N = 12.011, 1.008, 15.999, 14.007   # g/mol
    M_fuel = 3 * M_C + 8 * M_H
    M_air_per_O2 = 2 * M_O + 3.76 * 2 * M_N              # one "mole of O2" worth of air
    return M_fuel / (5 * M_air_per_O2)


def compressor_exit(T2, P2, PR, eta_c, gamma=1.4):
    """Temperature and pressure leaving the compressor.

    An ideal compressor raises T by the factor PR^((γ-1)/γ); a real one wastes
    work as friction, so the actual temperature rise is the ideal rise / η_c.
    """
    T3 = T2 * (1 + (PR ** ((gamma - 1) / gamma) - 1) / eta_c)
    P3 = PR * P2
    return T3, P3


def overall_far_hand(T3, T4, cp, LHV, eta_b):
    """Fuel-air ratio needed to heat the air from T3 to T4 (energy balance).

    Heat released by the fuel, f·LHV·η_b, must heat the whole gas (air + fuel,
    so 1+f) by cp·(T4−T3). Solving exactly for f gives f = cp·ΔT/(η_b·LHV − cp·ΔT).
    """
    q = cp * (T4 - T3)
    return q / (eta_b * LHV - q)


def equilibrate(phi, T, P):
    """Mix fuel and air at equivalence ratio φ, temperature T and pressure P, then
    let it react to chemical equilibrium at constant enthalpy and pressure (an
    adiabatic flame). Returns the gas object in its burned state."""
    gas.set_equivalence_ratio(phi, FUEL, OXIDIZER)  # sets the fuel/air mole fractions
    gas.TP = T, P                                  # sets the unburned-mixture state
    gas.equilibrate("HP")                          # burn: hold H and P constant
    return gas


def phi_overall_cantera(T3, P3, T4_target):
    """Equivalence ratio whose adiabatic equilibrium temperature equals T4_target.

    brentq searches between φ = 0.05 and 1.0 for the φ where T_eq(φ) − T4 = 0.
    Cantera assumes all fuel burns (100 % efficiency), so this comes out a bit
    LOWER than the hand value, which assumes only 95 % of the fuel's energy is released.
    """
    return brentq(lambda phi: equilibrate(phi, T3, P3).T - T4_target, 0.05, 1.0)


def air_split(phi_overall, phi_primary, phi_secondary, cooling_fraction):
    """How the total air is divided between the combustor zones.

    All the fuel enters the primary zone, so a zone that has received fraction X
    of the total air has φ_zone = φ_overall / X, i.e. X = φ_overall / φ_zone.
    Cooling air bypasses the zones along the liner wall and rejoins before the turbine.
    """
    X_pri = phi_overall / phi_primary            # cumulative air reaching primary
    X_sec = phi_overall / phi_secondary          # cumulative air reaching secondary
    split = {
        "primary": X_pri,
        "secondary": X_sec - X_pri,
        "dilution": 1.0 - X_sec - cooling_fraction,
        "cooling": cooling_fraction,
        "cum_primary": X_pri,
        "cum_secondary": X_sec,
    }
    for zone in ("primary", "secondary", "dilution"):
        assert split[zone] > 0, f"{zone} air fraction is not positive: {split[zone]:.4f}"
    assert split["cooling"] >= 0, "cooling_fraction cannot be negative"
    total = split["primary"] + split["secondary"] + split["dilution"] + split["cooling"]
    assert abs(total - 1.0) < 1e-9, f"air fractions sum to {total}, not 1"
    return split


def zone_temp_hand(phi, T3, cp, LHV, eta_b, far_st):
    """Zone temperature from a simple energy balance with an average cp.

    The fuel releases f·LHV·η_b per kg of air, which heats (1+f) kg of gas, so
    T = T3 + f·LHV·η_b / ((1+f)·cp). Ignores dissociation, so it overestimates hot zones.
    """
    assert phi <= 1.0, "hand formula assumes all fuel burns, valid only for φ ≤ 1"
    f = phi * far_st
    return T3 + f * LHV * eta_b / ((1 + f) * cp)


def zone_temp_equilibrium(phi, T3, P3):
    """Zone temperature and CO mole fraction from Cantera chemical equilibrium.

    Unlike the hand formula, equilibrium includes dissociation: at very high
    temperature some CO2 and H2O split back into CO, H2, OH etc., which absorbs
    heat and lowers the temperature.
    """
    # Why we can equilibrate the CUMULATIVE mixture from scratch instead of
    # mixing earlier-zone products with fresh air: the secondary zone's contents
    # are "primary-zone products + extra air". Mixing conserves total enthalpy
    # and total atoms, and an adiabatic equilibrium state depends only on the
    # atoms, enthalpy and pressure, not on the path taken. So "fuel + all the air
    # so far, burned at once" reaches the same equilibrium state.
    #
    # Note: in the cool dilution zone, real CO "freezes": reactions get too
    # slow to finish burning it, so real CO stays higher than the equilibrium
    # value. Treat the dilution-zone CO below as a lower bound.
    g = equilibrate(phi, T3, P3)
    return g.T, g["CO"].X[0]


def benchmark_check():
    """Check Cantera against a known result: stoichiometric propane-air at
    298.15 K and 1 atm has an adiabatic flame temperature of about 2260 K."""
    T_ad = equilibrate(1.0, 298.15, ct.one_atm).T
    ok = abs(T_ad - 2260.0) <= 30.0
    if ok:
        print(f"BENCHMARK: propane-air φ=1, 298.15 K, 1 atm -> {T_ad:.1f} K "
              f"(reference 2260 ± 30 K)  PASS\n")
    else:
        print("!" * 72)
        print(f"!!! BENCHMARK FAIL: got {T_ad:.1f} K, expected 2260 ± 30 K.")
        print("!!! Cantera / mechanism setup is wrong. Do NOT trust results below.")
        print("!" * 72 + "\n")
    return T_ad, ok


def sweep_and_plot(T3, P3, design_points, path):
    """Equilibrium temperature versus φ for the combustor inlet conditions.

    Shows how flame temperature peaks slightly rich of φ = 1 and falls on both
    sides; the marked points are the design zones. The sweep starts at φ = 0.1
    (not 0.4) so the curve also passes through the overall φ that sets T4.
    """
    phis = np.arange(0.1, 1.6 + 1e-9, 0.02)
    temps = np.array([equilibrate(p, T3, P3).T for p in phis])

    ink, muted, grid, line = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6"
    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
    ax.plot(phis, temps, color=line, linewidth=2)
    ax.axhline(T_turbine_max, color=muted, linewidth=1, linestyle="--")
    ax.text(1.6, T_turbine_max + 25, f"turbine limit {T_turbine_max:.0f} K",
            color=muted, fontsize=9, ha="right")
    for label, phi, T in design_points:
        ax.plot(phi, T, "o", markersize=8, color=line,
                markeredgecolor="white", markeredgewidth=2)
        ax.annotate(f"{label}\nφ={phi:.3g}, {T:.0f} K", (phi, T),
                    textcoords="offset points", xytext=(10, -4), fontsize=9,
                    color=ink, va="top")
    ax.set_xlabel("Equivalence ratio φ", color=ink)
    ax.set_ylabel("Equilibrium temperature (K)", color=ink)
    ax.set_title(f"Propane-air equilibrium temperature at T3 = {T3:.0f} K, "
                 f"P3 = {P3 / 1000:.0f} kPa", color=ink, fontsize=11)
    ax.grid(True, color=grid, linewidth=0.8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(muted)
    ax.tick_params(colors=muted)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def write_report(lines, path):
    """Print the report to the terminal and save the same text as Markdown."""
    text = "\n".join(lines)
    print(text)
    path.write_text(text + "\n", encoding="utf-8")


def md_table(header, rows):
    """Format a list of rows as a Markdown table (also readable in a terminal)."""
    widths = [max(len(str(x)) for x in col) for col in zip(header, *rows)]
    fmt = lambda r: "| " + " | ".join(str(x).ljust(w) for x, w in zip(r, widths)) + " |"
    out = [fmt(header), "|" + "|".join("-" * (w + 2) for w in widths) + "|"]
    out += [fmt(r) for r in rows]
    return out


def main():
    """Run every calculation and report the results."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    T_bench, bench_ok = benchmark_check()

    # --- Stations and fuel-air ratios ---------------------------------------
    far_st = stoich_far()
    T3, P3 = compressor_exit(T2, P2, PR, eta_c)

    f_hand = overall_far_hand(T3, T4_target, cp_overall, LHV, eta_b_exit)
    phi_hand = f_hand / far_st
    phi_cant = phi_overall_cantera(T3, P3, T4_target)
    f_cant = phi_cant * far_st

    if phi_overall_design is None:
        phi_used, phi_source = phi_hand, "hand-calculated (phi_overall_design = None)"
    else:
        phi_used, phi_source = phi_overall_design, "design input"
    f_used = phi_used * far_st
    m_fuel = f_used * m_air

    # T4 produced by the φ actually used, two ways
    T4_hand = zone_temp_hand(phi_used, T3, cp_overall, LHV, eta_b_exit, far_st)
    T4_cant, _ = zone_temp_equilibrium(phi_used, T3, P3)

    # --- Zones ----------------------------------------------------------------
    split = air_split(phi_used, phi_primary, phi_secondary, cooling_fraction)
    zones = [
        # name, φ, cumulative air, added air, cooling air (film along this zone's wall), cp
        ("Primary", phi_primary, split["cum_primary"], split["primary"],
         split["cooling"] * cooling_distribution["primary"], cp_primary),
        ("Secondary", phi_secondary, split["cum_secondary"], split["secondary"],
         split["cooling"] * cooling_distribution["secondary"], cp_secondary),
        ("Dilution", phi_used, 1.0, split["dilution"],
         split["cooling"] * cooling_distribution["dilution"], cp_dilution),
    ]
    zone_results = []
    for name, phi, cum, added, cool, cp in zones:
        T_h = zone_temp_hand(phi, T3, cp, LHV, eta_b_hot, far_st)
        T_c, X_co = zone_temp_equilibrium(phi, T3, P3)
        zone_results.append((name, phi, cum, added, cool, T_h, T_c, X_co))

    # --- Sanity checks ----------------------------------------------------------
    checks = []
    _, _, _, _, _, Th_p, Tc_p, _ = zone_results[0]
    if Tc_p < Th_p:
        checks.append(f"PASS  Primary: Cantera {Tc_p:.0f} K < hand {Th_p:.0f} K "
                      f"(dissociation lowers it by {Th_p - Tc_p:.0f} K)")
    else:
        checks.append(f"WARN  Primary: Cantera {Tc_p:.0f} K is NOT below hand {Th_p:.0f} K "
                      f"(expected dissociation; check cp_primary)")

    _, _, _, _, _, Th_s, Tc_s, _ = zone_results[1]
    d_s = Th_s - Tc_s
    if abs(d_s) < 30:
        checks.append(f"INFO  Secondary: hand - Cantera = {d_s:+.0f} K (agree within 30 K; "
                      f"at lean φ dissociation is small, so this mainly tests cp_secondary)")
    elif d_s > 0:
        checks.append(f"PASS  Secondary: Cantera {Tc_s:.0f} K < hand {Th_s:.0f} K "
                      f"({d_s:.0f} K lower)")
    else:
        checks.append(f"WARN  Secondary: hand is {-d_s:.0f} K BELOW Cantera "
                      f"(cp_secondary is probably too high)")

    air_total = split["primary"] + split["secondary"] + split["dilution"] + split["cooling"]
    checks.append(f"{'PASS' if abs(air_total - 1) < 1e-9 else 'WARN'}  "
                  f"Air fractions sum to {air_total * 100:.2f} %")

    for label, T4 in (("hand (η_b = %.2f)" % eta_b_exit, T4_hand),
                      ("Cantera (ideal, upper bound)", T4_cant)):
        margin = 1 - T4 / T_turbine_max
        if T4 >= T_turbine_max:
            status = "FAIL"
        elif margin < T4_margin_min:
            status = "WARN"
        else:
            status = "PASS"
        checks.append(f"{status}  T4 {label}: {T4:.1f} K ({c_of(T4):.1f} °C), "
                      f"margin {margin * 100:.1f} % below {T_turbine_max:.0f} K limit "
                      f"(need ≥ {T4_margin_min * 100:.0f} %)")

    # --- Plot -------------------------------------------------------------------
    plot_path = OUT_DIR / "T_vs_phi.png"
    sweep_and_plot(T3, P3, [
        ("Primary", phi_primary, zone_results[0][6]),
        ("Secondary", phi_secondary, zone_results[1][6]),
        ("Overall (T4)", phi_used, T4_cant),
    ], plot_path)

    # --- Report -------------------------------------------------------------------
    L = []
    L.append(f"# Combustion chamber calculation - {date.today().isoformat()}")
    L.append("")
    L.append(f"Benchmark (propane-air φ=1, 298.15 K, 1 atm): {T_bench:.1f} K "
             f"-> {'PASS' if bench_ok else 'FAIL'}")
    L.append("")
    L.append("## Inputs")
    L.append("")
    L += md_table(["Input", "Value"], [
        ["m_air", f"{m_air} kg/s"],
        ["T2 / P2", f"{T2} K / {P2} Pa"],
        ["PR / eta_c", f"{PR} / {eta_c}"],
        ["T4_target", f"{T4_target} K ({c_of(T4_target):.1f} °C)"],
        ["T_turbine_max / margin", f"{T_turbine_max} K / {T4_margin_min * 100:.0f} %"],
        ["phi_overall_design", f"{phi_overall_design}"],
        ["phi_primary / phi_secondary", f"{phi_primary} / {phi_secondary}"],
        ["cooling_fraction", f"{cooling_fraction}"],
        ["LHV", f"{LHV / 1e6:.2f} MJ/kg"],
        ["eta_b_exit / eta_b_hot", f"{eta_b_exit} / {eta_b_hot}"],
        ["cp overall / pri / sec / dil", f"{cp_overall} / {cp_primary} / {cp_secondary} / {cp_dilution} J/kg·K"],
        ["Mechanism / fuel", f"{MECH} / {FUEL}"],
    ])
    L.append("")
    L.append("## Stations")
    L.append("")
    L += md_table(["Quantity", "Value"], [
        ["FAR_st", f"{far_st:.5f}"],
        ["T3", f"{T3:.2f} K ({c_of(T3):.2f} °C)"],
        ["P3", f"{P3:,.0f} Pa ({P3 / 1e5:.3f} bar)"],
    ])
    L.append("")
    L.append("## Overall fuel-air ratio (to reach T4_target)")
    L.append("")
    L += md_table(["Method", "f", "φ_overall"], [
        [f"Hand (cp={cp_overall}, η_b={eta_b_exit})", f"{f_hand:.5f}", f"{phi_hand:.4f}"],
        ["Cantera (equilibrium, η_b=1)", f"{f_cant:.5f}", f"{phi_cant:.4f}"],
    ])
    L.append("")
    L.append(f"Cantera needs {(1 - phi_cant / phi_hand) * 100:.1f} % less fuel: it assumes every "
             f"bit of fuel burns (η_b = 1) and uses exact cp(T), while the hand method "
             f"assumes 95 % efficiency and an average cp of {cp_overall} J/kg·K.")
    L.append("")
    L.append("## Design point")
    L.append("")
    L += md_table(["Quantity", "Value"], [
        ["φ_overall used", f"{phi_used:.4f} ({phi_source})"],
        ["f used", f"{f_used:.5f}"],
        ["Fuel mass flow", f"{m_fuel * 1000:.3f} g/s"],
        [f"T4 hand (η_b={eta_b_exit})", f"{T4_hand:.1f} K ({c_of(T4_hand):.1f} °C)"],
        ["T4 Cantera (ideal)", f"{T4_cant:.1f} K ({c_of(T4_cant):.1f} °C)"],
    ])
    L.append("")
    L.append("## Zones")
    L.append("")
    rows = []
    for name, phi, cum, added, cool, T_h, T_c, X_co in zone_results:
        rows.append([
            name, f"{phi:.3f}", f"{cum * 100:.1f}", f"{added * 100:.1f}", f"{cool * 100:.1f}",
            f"{T_h:.0f} K ({c_of(T_h):.0f} °C)", f"{T_c:.0f} K ({c_of(T_c):.0f} °C)",
            f"{T_h - T_c:+.0f}", f"{X_co * 100:.3g} % ({X_co * 1e6:.3g} ppm)",
        ])
    L += md_table(["Zone", "φ", "Cum. air %", "Air added %", "Cooling %",
                   f"T hand (η_b={eta_b_hot})", "T Cantera", "Hand−Cantera (K)", "CO (equilibrium)"],
                  rows)
    L.append("")
    L.append("Air added % = air entering through that zone's holes. Cooling % = wall-film air "
             "running along that zone's liner wall (% of total air); it is NOT part of the "
             "zone's φ and mixes in only before the turbine.")
    L.append("")
    L.append("Dilution-zone CO is an equilibrium lower bound; real CO freezes higher.")
    L.append("")
    L.append("## Sanity checks")
    L.append("")
    L += [f"- {c}" for c in checks]
    L.append("")
    L.append(f"Plot saved: {plot_path.relative_to(Path(__file__).parent.parent).as_posix()}")
    L += liner_section(T3, P3, split, phi_used, far_st,  # added: liner dP, hole and cooling areas
                       zone_results, T4_hand, T4_cant)

    write_report(L, OUT_DIR / "chamber_calc_results.md")


# =============================================================================
# Liner pressure drop and hole areas  (added below the existing functions)
# =============================================================================
HOLE_ZONES = ("primary", "secondary", "dilution", "cooling")


def air_density(P, T, R):
    """Air density from the ideal-gas law, ρ = P/(R·T).

    Evaluated at T3 and P3: the air sitting in the annulus just before it
    passes through the liner holes.
    """
    return P / (R * T)


def jet_velocity(rho, dP):
    """Ideal speed of the air jet through a liner hole, V = sqrt(2·ΔP/ρ).

    Bernoulli: the whole liner pressure drop is turned into jet kinetic energy,
    ΔP = ½·ρ·V². Higher ΔP means faster jets that penetrate deeper into the hot gas.
    """
    return np.sqrt(2 * dP / rho)


def hole_area_for_flow(mdot, rho, dP, Cd):
    """Geometric hole area needed to pass mass flow mdot, A = mdot / (Cd·sqrt(2·ρ·ΔP)).

    From mdot = Cd·A·sqrt(2·ρ·ΔP). Cd corrects for the jet narrowing after it
    leaves the hole (vena contracta), so the real hole must be larger than the
    ideal flow area.
    """
    return mdot / (Cd * np.sqrt(2 * rho * dP))


def zone_hole_areas(m_air, split, rho, dP, Cd):
    """Hole area per zone (m²), each sized for the air ADDED in that zone.

    Every hole sees the same annulus air (same ρ) and the same liner ΔP, so the
    areas are in the same proportion as the air shares from air_split().
    Cooling is included only when cooling_fraction > 0.
    """
    areas = {}
    for zone in HOLE_ZONES:
        if zone == "cooling" and not split["cooling"] > 0:
            continue
        areas[zone] = hole_area_for_flow(split[zone] * m_air, rho, dP, Cd)
    return areas


def candidate_table(A_zone, d_lo=None, d_hi=None, n_max=24):
    """Reference table of possible hole patterns that give area A_zone (m²).

    For each even hole count n from 4 to n_max: the exact diameter for n equal round
    holes, d = sqrt(4·A/(n·π)), the nearest drillable diameter on the d_step_mm
    grid, and the area error that rounding causes. Status is OK when the drilled
    diameter lies inside [d_lo, d_hi] (default [d_min_mm, d_max_mm]), else OUT. If
    liner_diameter_mm is set, also gives how many holes fit around one circumference
    at the minimum spacing and how many rows n holes would need. Does NOT choose a pattern.
    """
    d_lo = d_min_mm if d_lo is None else d_lo
    d_hi = d_max_mm if d_hi is None else d_hi
    rows = []
    for n in range(4, n_max + 1, 2):
        d_exact = np.sqrt(4 * A_zone / (n * np.pi)) * 1000          # mm
        d_drill = np.floor(d_exact / d_step_mm + 0.5) * d_step_mm   # nearest grid value, mm
        A_drill = n * np.pi / 4 * (d_drill / 1000) ** 2
        err = (A_drill / A_zone - 1) * 100
        status = "OK" if d_lo <= d_drill <= d_hi else "OUT"
        row = {"n": n, "d_exact": d_exact, "d_drill": d_drill, "err": err, "status": status}
        if liner_diameter_mm is not None:
            per_row = int(np.floor(np.pi * liner_diameter_mm / (spacing_factor * d_drill)))
            row["per_row"] = per_row
            row["rows"] = int(np.ceil(n / per_row)) if per_row > 0 else None
        rows.append(row)
    return rows


def cooling_hole_areas(m_air, split, rho, dP, Cd_main, Cd_alt):
    """Cooling-hole area (m²) for the film air along each zone's wall.

    The cooling air (cooling_fraction × m_air) is divided by cooling_distribution
    and sized with the same ρ and ΔP as the zone holes, A = mdot/(Cd·sqrt(2·ρ·ΔP)).
    Returns {group: (area at Cd_main, area at Cd_alt)}; a higher Cd needs less area.
    """
    return {z: (hole_area_for_flow(split["cooling"] * share * m_air, rho, dP, Cd_main),
                hole_area_for_flow(split["cooling"] * share * m_air, rho, dP, Cd_alt))
            for z, share in cooling_distribution.items()}


def liner_section(T3, P3, split, phi_used, far_st, zone_results, T4_hand, T4_cant):
    """Report lines for the liner pressure drop, hole-area and wall-cooling section,
    including its PASS/WARN/INFO checks. Uses results already computed in main()."""
    rho3 = air_density(P3, T3, R_air)
    cm2 = 1e4  # m² -> cm²
    cf = split["cooling"]
    zones = ("primary", "secondary", "dilution")   # jet holes, sized with Cd
    groups = zones + (("cooling",) if cf > 0 else ())

    def areas_at(dP, sp):
        """Zone-hole areas (Cd) plus total cooling-hole area (Cd_cooling) at one dP."""
        a = zone_hole_areas(m_air, sp, rho3, dP, Cd)
        cool = cooling_hole_areas(m_air, sp, rho3, dP, Cd_cooling, Cd_cooling_sensitivity)
        a["cooling"] = sum(v[0] for v in cool.values()) if sp["cooling"] > 0 else 0.0
        return a, cool

    def total_area(dP, cf_):
        """Total hole area computed independently: all zone-hole air at Cd + cooling air at Cd_cooling."""
        return (hole_area_for_flow((1 - cf_) * m_air, rho3, dP, Cd)
                + (hole_area_for_flow(cf_ * m_air, rho3, dP, Cd_cooling) if cf_ > 0 else 0.0))

    L = [""]
    L.append("## Liner pressure drop and hole areas")
    L.append("")
    L += md_table(["Input", "Value"], [
        ["P3", f"{P3:,.0f} Pa"],
        ["T3", f"{T3:.2f} K"],
        ["rho3 = P3/(R_air·T3)", f"{rho3:.3f} kg/m³ (R_air = {R_air} J/kg·K)"],
        ["Cd", f"{Cd} (assumed, sharp-edged hole)"],
        ["Design dP fraction", f"{dP_liner_fraction * 100:g} % of P3"],
        ["cooling_fraction", f"{cf * 100:g} % of total air"],
        ["cooling_distribution (provisional)",
         " / ".join(f"{z} {s * 100:g} %" for z, s in cooling_distribution.items())],
        ["Cd_cooling / sensitivity", f"{Cd_cooling} (placeholder) / {Cd_cooling_sensitivity}"],
        ["Cooling hole limits", f"{d_cool_min_mm} - {d_cool_max_mm} mm"],
    ])
    L.append("")
    L.append("Assumption: cooling air runs along the inner liner wall as a film and does NOT take "
             "part in any zone's φ; it mixes into the gas before the turbine. Primary/secondary φ, "
             "cumulative air and temperatures are therefore unchanged; only the dilution-hole air "
             "shrinks. φ_overall and T4 are unchanged (same total air and fuel).")
    L.append("")

    # --- dP comparison -------------------------------------------------------
    L.append("### dP comparison")
    L.append("")
    rows = []
    for frac in dP_sweep:
        dP = frac * P3
        areas, _ = areas_at(dP, split)
        A_tot = total_area(dP, cf)
        mark = " (design)" if abs(frac - dP_liner_fraction) < 1e-12 else ""
        rows.append([f"{frac * 100:g} %{mark}", f"{dP:,.0f}", f"{dP / 1000:.2f}",
                     f"{jet_velocity(rho3, dP):.1f}", f"{A_tot * cm2:.2f}"]
                    + [f"{areas[z] * cm2:.2f}" for z in groups])
    L += md_table(["dP / P3", "dP (Pa)", "dP (kPa)", "Jet velocity (m/s)", "Total area (cm²)"]
                  + [f"{z.capitalize()} (cm²)" for z in groups], rows)
    L.append("")

    # --- Design point --------------------------------------------------------
    dP = dP_liner_fraction * P3
    areas, cool = areas_at(dP, split)
    A_tot = total_area(dP, cf)
    L.append(f"### Design point (dP = {dP_liner_fraction * 100:g} % of P3 = {dP:,.0f} Pa, "
             f"jet velocity {jet_velocity(rho3, dP):.1f} m/s)")
    L.append("")
    labels = {"primary": "Primary holes", "secondary": "Secondary holes",
              "dilution": "Dilution holes", "cooling": f"Cooling, all walls (Cd {Cd_cooling})"}
    rows = [[labels[z], f"{split[z] * 100:.1f}", f"{split[z] * m_air:.4f}",
             f"{areas[z] * cm2:.2f}"] for z in groups]
    share_sum = sum(split[z] for z in groups)
    A_sum = sum(areas[z] for z in groups)
    rows.append(["Total", f"{share_sum * 100:.1f}", f"{share_sum * m_air:.4f}", f"{A_sum * cm2:.2f}"])
    L += md_table(["Group", "Air share %", "Air flow (kg/s)", "Hole area (cm²)"], rows)
    L.append("")

    # --- Cooling holes -------------------------------------------------------
    if cf > 0:
        L.append("### Wall cooling holes at design dP")
        L.append("")
        rows = [[f"{z.capitalize()} wall", f"{s * 100:g}", f"{cf * s * 100:.1f}",
                 f"{cf * s * m_air:.4f}", f"{cool[z][0] * cm2:.2f}", f"{cool[z][1] * cm2:.2f}"]
                for z, s in cooling_distribution.items()]
        rows.append(["Total", f"{sum(cooling_distribution.values()) * 100:g}", f"{cf * 100:.1f}",
                     f"{cf * m_air:.4f}", f"{sum(v[0] for v in cool.values()) * cm2:.2f}",
                     f"{sum(v[1] for v in cool.values()) * cm2:.2f}"])
        L += md_table(["Group", "% of cooling air", "% of total air", "Air flow (kg/s)",
                       f"Area at Cd {Cd_cooling} (cm²)",
                       f"Area at Cd {Cd_cooling_sensitivity} (cm²)"], rows)
        L.append("")

    # --- Cooling 0 % vs current ----------------------------------------------
    L.append(f"### Cooling 0 % vs {cf * 100:g} % (design dP)")
    L.append("")
    rows = []
    for c in sorted({0.0, cf}):
        sp = air_split(phi_used, phi_primary, phi_secondary, c)
        a, _ = areas_at(dP, sp)
        T4_h = zone_temp_hand(phi_used, T3, cp_overall, LHV, eta_b_exit, far_st)
        rows.append([f"{c * 100:g} %", f"{sp['dilution'] * 100:.1f}", f"{a['dilution'] * cm2:.2f}",
                     f"{a['cooling'] * cm2:.2f}", f"{total_area(dP, c) * cm2:.2f}",
                     f"{phi_used:.4f}", f"{T4_h:.1f}"])
    L += md_table(["Cooling", "Dilution-hole air %", "Dilution hole area (cm²)",
                   "Cooling area (cm²)", "Total area (cm²)", "φ_overall", "T4 hand (K)"], rows)
    L.append("")
    L.append("Plain radial holes produce jets, not a wall film. Cooling holes need a deflector lip, "
             "splash ring, or angled drilling; geometry not designed yet.")
    L.append("")

    # --- Candidate tables ----------------------------------------------------
    L.append("### Hole candidates at design dP (reference only, nothing chosen)")
    L.append("")
    L.append(f"Drill grid {d_step_mm} mm; OK = drilled diameter within "
             f"[{d_min_mm}, {d_max_mm}] mm. Area error = drilled area vs required area.")
    if liner_diameter_mm is None:
        L.append("Per-row fit skipped: liner_diameter_mm = None (placeholder, not chosen yet).")
    else:
        L.append(f"Per-row fit: liner diameter {liner_diameter_mm} mm, "
                 f"min spacing {spacing_factor} × d.")

    def add_candidates(title, A, **limits):
        L.append("")
        L.append(f"#### {title} (required area {A * cm2:.2f} cm²)")
        L.append("")
        header = ["Holes", "d exact (mm)", "d drilled (mm)", "Area error", "Status"]
        if liner_diameter_mm is not None:
            header += ["Max holes/row", "Rows needed"]
        rows = []
        for r in candidate_table(A, **limits):
            row = [r["n"], f"{r['d_exact']:.2f}", f"{r['d_drill']:.1f}",
                   f"{r['err']:+.1f} %", r["status"]]
            if liner_diameter_mm is not None:
                row += [r["per_row"], "—" if r["rows"] is None else r["rows"]]
            rows.append(row)
        L.extend(md_table(header, rows))

    for z in zones:
        add_candidates(z.capitalize(), areas[z])
    L.append("")

    if cf > 0:
        L.append(f"### Cooling hole candidates at design dP, Cd {Cd_cooling} "
                 f"(reference only, nothing chosen)")
        L.append("")
        L.append(f"Even counts 4-60; OK = drilled diameter within "
                 f"[{d_cool_min_mm}, {d_cool_max_mm}] mm.")
        for z in cooling_distribution:
            add_candidates(f"{z.capitalize()} wall cooling", cool[z][0],
                           d_lo=d_cool_min_mm, d_hi=d_cool_max_mm, n_max=60)
        L.append("")

    # --- Checks --------------------------------------------------------------
    checks = []
    sum_err = abs(A_sum / A_tot - 1) * 100
    checks.append(f"{'PASS' if sum_err < 0.1 else 'WARN'}  Group areas sum to "
                  f"{A_sum * cm2:.3f} cm² vs total {A_tot * cm2:.3f} cm² "
                  f"(difference {sum_err:.2e} %, limit 0.1 %)")
    for z in groups:
        Cd_z = Cd_cooling if z == "cooling" else Cd
        area_ratio = areas[z] / hole_area_for_flow(m_air, rho3, dP, Cd_z)
        ratio_err = abs(area_ratio / split[z] - 1) * 100
        checks.append(f"{'PASS' if ratio_err < 0.1 else 'WARN'}  {z.capitalize()}: area share "
                      f"{area_ratio * 100:.2f} % vs air share {split[z] * 100:.2f} % at Cd {Cd_z} "
                      f"(difference {ratio_err:.2e} %, limit 0.1 %)")

    dist_sum = sum(cooling_distribution.values())
    checks.append(f"{'PASS' if abs(dist_sum - 1) < 1e-9 else 'WARN'}  cooling_distribution "
                  f"sums to {dist_sum * 100:.2f} % (must be 100 %)")
    checks.append(f"{'PASS' if split['dilution'] > 0 else 'WARN'}  Dilution-hole air share "
                  f"{split['dilution'] * 100:.1f} % > 0 (also asserted in air_split)")
    shares = [split[z] for z in zones] + [cf * s for s in cooling_distribution.values()]
    checks.append(f"{'PASS' if abs(sum(shares) - 1) < 1e-9 else 'WARN'}  Zone-hole + cooling "
                  f"air shares sum to {sum(shares) * 100:.2f} %")

    # Same quantities recomputed with cooling_fraction = 0
    sp0 = air_split(phi_used, phi_primary, phi_secondary, 0.0)
    now, ref = [], []
    for i, z in enumerate(("primary", "secondary")):
        _, phi_z, cum, added, _, T_h, T_c, _ = zone_results[i]
        T_h0 = zone_temp_hand((phi_primary, phi_secondary)[i], T3, (cp_primary, cp_secondary)[i],
                              LHV, eta_b_hot, far_st)
        T_c0, _ = zone_temp_equilibrium((phi_primary, phi_secondary)[i], T3, P3)
        now += [phi_z, cum, added, T_h, T_c]
        ref += [(phi_primary, phi_secondary)[i], sp0[f"cum_{z}"], sp0[z], T_h0, T_c0]
    T4_hand_0 = zone_temp_hand(phi_used, T3, cp_overall, LHV, eta_b_exit, far_st)
    T4_cant_0, _ = zone_temp_equilibrium(phi_used, T3, P3)
    now += [T4_hand, T4_cant]
    ref += [T4_hand_0, T4_cant_0]
    same = all(abs(a - b) <= 1e-9 * max(1.0, abs(b)) for a, b in zip(now, ref))
    checks.append(f"{'PASS' if same else 'WARN'}  φ_overall {phi_used:.4f}, T4 hand {T4_hand:.1f} K, "
                  f"T4 Cantera {T4_cant:.1f} K, primary/secondary φ, cumulative air % and "
                  f"temperatures identical to the 0 % cooling run")

    A_same_Cd = (sum(areas[z] for z in zones)
                 + (hole_area_for_flow(cf * m_air, rho3, dP, Cd) if cf > 0 else 0.0))
    A_0 = total_area(dP, 0.0)
    tot_err = abs(A_same_Cd / A_0 - 1) * 100
    checks.append(f"{'PASS' if tot_err < 0.1 else 'WARN'}  Total hole area at the same Cd ({Cd}): "
                  f"{A_same_Cd * cm2:.2f} cm² vs {A_0 * cm2:.2f} cm² with 0 % cooling "
                  f"(difference {tot_err:.2e} %, limit 0.1 %)")

    if liner_diameter_mm is None:
        checks.append("INFO  Hole counts and diameters not chosen: liner diameter not set.")
    else:
        checks.append("INFO  Hole counts and diameters not chosen: candidate tables are reference only.")
    checks.append("INFO  Cooling hole geometry (film-forming lip / splash ring / angled holes) "
                  "not designed.")
    L.append("### Liner checks")
    L.append("")
    L += [f"- {c}" for c in checks]
    return L


if __name__ == "__main__":
    main()
