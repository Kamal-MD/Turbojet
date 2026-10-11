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
cooling_fraction = 0.0      # -, fraction of total air for liner wall cooling (not yet allocated)
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
        # name, φ, cumulative air, added air, cooling air, cp
        ("Primary", phi_primary, split["cum_primary"], split["primary"], 0.0, cp_primary),
        ("Secondary", phi_secondary, split["cum_secondary"], split["secondary"], 0.0, cp_secondary),
        ("Dilution", phi_used, 1.0, split["dilution"], split["cooling"], cp_dilution),
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
    L.append("Dilution-zone CO is an equilibrium lower bound; real CO freezes higher.")
    L.append("")
    L.append("## Sanity checks")
    L.append("")
    L += [f"- {c}" for c in checks]
    L.append("")
    L.append(f"Plot saved: {plot_path.relative_to(Path(__file__).parent.parent).as_posix()}")

    write_report(L, OUT_DIR / "chamber_calc_results.md")


if __name__ == "__main__":
    main()
