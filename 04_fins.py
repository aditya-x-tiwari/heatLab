import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# USER-EDITABLE CONSTANTS (from the lab sheet)
# ============================================================
DATA_FILE = "input_data.csv"

D_FIN_MM = 12.7          # fin diameter in mm
L_FIN_MM = 150.0         # fin length in mm
D_ORIFICE_MM = 10.0      # orifice diameter in mm
DUCT_W_MM = 100.0        # duct width in mm
DUCT_H_MM = 150.0        # duct height in mm
CD_ORIFICE = 0.62        # coefficient of discharge
K_FIN_W_MK = 111.0       # brass fin thermal conductivity (editable)
P_ATM_PA = 101325.0      # assumed atmospheric pressure for air-property evaluation
R_AIR = 287.05           # J/(kg.K)
CP_AIR = 1007.0          # J/(kg.K)
PR_AIR = 0.71            # approx. for air
RHO_WATER = 1000.0       # kg/m^3
G = 9.81                 # m/s^2

# Thermocouple positions along the fin.
# The sheet does not list exact locations, so equal spacing is used by default.
# Edit these if your apparatus uses different TC positions.
X_POS_MM = np.linspace(0.0, L_FIN_MM, 5)


# ============================================================
# AIR PROPERTIES
# ============================================================
def air_properties(temp_c: float, pressure_pa: float = P_ATM_PA):
    """
    Return air properties at the given temperature.
    Uses Sutherland's law for viscosity and ideal-gas density.
    Thermal conductivity is estimated from mu*cp/Pr.
    """
    t_k = temp_c + 273.15

    # Sutherland's law for dynamic viscosity
    mu0 = 1.716e-5   # Pa.s at T0
    t0 = 273.15      # K
    s = 111.0        # K
    mu = mu0 * ((t_k / t0) ** 1.5) * ((t0 + s) / (t_k + s))

    rho = pressure_pa / (R_AIR * t_k)
    nu = mu / rho
    k_air = mu * CP_AIR / PR_AIR  # consistent with PR_AIR
    return {
        "rho": rho,
        "mu": mu,
        "nu": nu,
        "k": k_air,
        "pr": PR_AIR,
        "t_k": t_k,
    }


# ============================================================
# CORE CALCULATIONS
# ============================================================
def orifice_flow_rate(head_cm: float, ambient_temp_c: float) -> float:
    """
    Volumetric flow rate through the orifice meter, following the lab-sheet form:
        Q = Cd * (pi/4) * d_orifice^2 * sqrt(2 g h_w * rho_w / rho_a)

    head_cm is the manometer differential head in cm of water.
    """
    head_m = abs(head_cm) / 100.0
    rho_air = air_properties(ambient_temp_c)["rho"]

    d_orifice_m = D_ORIFICE_MM / 1000.0
    area_orifice = math.pi * (d_orifice_m ** 2) / 4.0

    q = CD_ORIFICE * area_orifice * math.sqrt(2.0 * G * head_m * (RHO_WATER / rho_air))
    return q


def velocity_at_mean_film(q_m3_s: float, ambient_temp_c: float, mean_film_temp_c: float) -> tuple[float, float]:
    """
    Returns:
      V_f  = duct air velocity at ambient temperature
      V_mf = air velocity corrected to mean film temperature
    """
    duct_area_m2 = (DUCT_W_MM / 1000.0) * (DUCT_H_MM / 1000.0)
    v_f = q_m3_s / duct_area_m2
    v_mf = v_f * ((mean_film_temp_c + 273.15) / (ambient_temp_c + 273.15))
    return v_f, v_mf


def nusselt_from_re(re: float) -> tuple[float, str]:
    """
    Correlation from the lab sheet for cross-flow over a cylinder / fin.
    """
    if 40.0 <= re < 4000.0:
        return 0.615 * (re ** 0.466), "Nu = 0.615 Re^0.466"
    elif 4000.0 <= re < 40000.0:
        return 0.174 * (re ** 0.618), "Nu = 0.174 Re^0.618"
    else:
        # Keep the code running even if outside the stated range.
        # Use the closer available fit by Reynolds range.
        if re < 40.0:
            return 0.615 * (re ** 0.466), "outside range; used low-Re fit"
        return 0.174 * (re ** 0.618), "outside range; used mid-Re fit"


def fin_temperature_profiles(
    x_m: np.ndarray,
    tb_c: float,
    t_inf_c: float,
    h: float,
    k_fin: float,
    d_fin_m: float,
    l_fin_m: float,
):
    """
    Returns temperature distributions for:
      1) insulated / adiabatic tip
      2) ideal tip (tip held at ambient temperature)
      3) convective tip

    Temperatures are returned in degC and the formulas are in terms of theta = T - T_inf.
    """
    p = math.pi * d_fin_m
    a_c = math.pi * (d_fin_m ** 2) / 4.0
    theta_b = tb_c - t_inf_c

    m = math.sqrt(h * p / (k_fin * a_c))
    m_l = m * l_fin_m
    beta = h / (m * k_fin)  # A_t/A_c = 1 for circular pin fin

    # Temperature distributions (theta/theta_b)
    # Insulated tip
    theta_ratio_ins = np.cosh(m * (l_fin_m - x_m)) / np.cosh(m_l)

    # Ideal tip (tip temperature equals ambient)
    # theta(L) = 0
    theta_ratio_ideal = np.exp( -m * x_m)

    # Convective tip
    theta_ratio_conv = (
        np.cosh(m * (l_fin_m - x_m)) + beta * np.sinh(m * (l_fin_m - x_m))
    ) / (
        np.cosh(m_l) + beta * np.sinh(m_l)
    )

    t_ins = t_inf_c + theta_b * theta_ratio_ins
    t_ideal = t_inf_c + theta_b * theta_ratio_ideal
    t_conv = t_inf_c + theta_b * theta_ratio_conv

    # Heat transfer rates
    q_ins = math.sqrt(h * p * k_fin * a_c) * theta_b * math.tanh(m_l)
    q_ideal = math.sqrt(h * p * k_fin * a_c) * theta_b / math.tanh(m_l)  # coth(mL)
    q_conv = math.sqrt(h * p * k_fin * a_c) * theta_b * (
        (math.sinh(m_l) + beta * math.cosh(m_l)) / (math.cosh(m_l) + beta * math.sinh(m_l))
    )

    # Fin efficiency (standard adiabatic-tip efficiency)
    eta_ins = math.tanh(m_l) / m_l if m_l != 0 else float("nan")

    # Fin effectiveness = q_actual / (h * A_c * theta_b)
    # This is a practical comparison for a pin fin.
    effectiveness_ins = q_ins / (h * a_c * theta_b) if theta_b != 0 else float("nan")
    effectiveness_ideal = q_ideal / (h * a_c * theta_b) if theta_b != 0 else float("nan")
    effectiveness_conv = q_conv / (h * a_c * theta_b) if theta_b != 0 else float("nan")

    return {
        "m": m,
        "beta": beta,
        "A_c": a_c,
        "P": p,
        "theta_b": theta_b,
        "t_ins": t_ins,
        "t_ideal": t_ideal,
        "t_conv": t_conv,
        "q_ins_W": q_ins,
        "q_ideal_W": q_ideal,
        "q_conv_W": q_conv,
        "eta_ins": eta_ins,
        "effectiveness_ins": effectiveness_ins,
        "effectiveness_ideal": effectiveness_ideal,
        "effectiveness_conv": effectiveness_conv,
    }


# ============================================================
# MAIN
# ============================================================
def main():
    data_path = Path(DATA_FILE)
    if not data_path.exists():
        raise FileNotFoundError(f"Could not find {DATA_FILE}. Put the CSV in the same folder as the script.")

    df = pd.read_csv(data_path)

    # Constants in SI
    d_fin_m = D_FIN_MM / 1000.0
    l_fin_m = L_FIN_MM / 1000.0
    x_m = X_POS_MM / 1000.0

    rows = []

    for _, row in df.iterrows():
        run = int(row["run"])
        v = float(row["voltage_V"])
        h1 = float(row["manometer_top_cm"])
        h2 = float(row["manometer_bottom_cm"])
        head_cm = float(row["manometer_diff_cm"])

        temps = np.array([row["T1_C"], row["T2_C"], row["T3_C"], row["T4_C"], row["T5_C"]], dtype=float)
        t_amb = float(row["ambient_C"])

        t_mean_fin = float(np.mean(temps))
        t_film = 0.5 * (t_mean_fin + t_amb)

        # Orifice flow + velocities
        qv = orifice_flow_rate(head_cm=head_cm, ambient_temp_c=t_amb)
        v_f, v_mf = velocity_at_mean_film(qv, ambient_temp_c=t_amb, mean_film_temp_c=t_film)

        # Air properties at mean film temperature
        air = air_properties(t_film)

        # Reynolds number and Nusselt number
        re = (air["rho"] * v_mf * d_fin_m) / air["mu"]
        nu, nu_fit_used = nusselt_from_re(re)
        h_conv = nu * air["k"] / d_fin_m

        # Fin temperatures and heat rates
        fin = fin_temperature_profiles(
            x_m=x_m,
            tb_c=float(temps[0]),
            t_inf_c=t_amb,
            h=h_conv,
            k_fin=K_FIN_W_MK,
            d_fin_m=d_fin_m,
            l_fin_m=l_fin_m,
        )

        rows.append({
            "run": run,
            "voltage_V": v,
            "manometer_top_cm": h1,
            "manometer_bottom_cm": h2,
            "manometer_diff_cm": head_cm,
            "ambient_C": t_amb,
            "T1_C": temps[0],
            "T2_C": temps[1],
            "T3_C": temps[2],
            "T4_C": temps[3],
            "T5_C": temps[4],
            "T_mean_fin_C": t_mean_fin,
            "T_film_C": t_film,
            "Qv_m3_s": qv,
            "V_duct_m_s": v_f,
            "V_mean_film_m_s": v_mf,
            "rho_air_kg_m3": air["rho"],
            "mu_air_Pa_s": air["mu"],
            "nu_air_m2_s": air["nu"],
            "k_air_W_mK": air["k"],
            "Re": re,
            "Nu": nu,
            "Nu_formula": nu_fit_used,
            "h_W_m2K": h_conv,
            "m_1_m": fin["m"],
            "beta": fin["beta"],
            "q_ins_W": fin["q_ins_W"],
            "q_ideal_W": fin["q_ideal_W"],
            "q_conv_W": fin["q_conv_W"],
            "eta_ins": fin["eta_ins"],
            "effectiveness_ins": fin["effectiveness_ins"],
            "effectiveness_ideal": fin["effectiveness_ideal"],
            "effectiveness_conv": fin["effectiveness_conv"],
        })

    results = pd.DataFrame(rows)
    results.to_csv("output_results.csv", index=False)

    # ============================================================
    # PRINT RESULTS
    # ============================================================
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 160)
    print("\n=== Computed Results ===\n")
    print(results[[
        "run", "manometer_diff_cm", "T_mean_fin_C", "T_film_C", "V_mean_film_m_s",
        "Re", "Nu", "h_W_m2K", "m_1_m", "q_ins_W", "q_ideal_W", "q_conv_W",
        "eta_ins", "effectiveness_ins", "effectiveness_ideal", "effectiveness_conv"
    ]].to_string(index=False))

    print("\nResults saved to output_results.csv")

    out_dir = Path("plots")
    out_dir.mkdir(exist_ok=True)

    # ============================================================
    # PLOT 1: Observed temperature distribution for each run
    # ============================================================
    fig, axes = plt.subplots(len(results), 1, figsize=(8, 12), sharex=True)
    if len(results) == 1:
        axes = [axes]

    for ax, (_, r) in zip(axes, results.iterrows()):
        obs_t = np.array([r["T1_C"], r["T2_C"], r["T3_C"], r["T4_C"], r["T5_C"]], dtype=float)
        ax.plot(X_POS_MM, obs_t, marker="o", label="Observed")
        ax.set_title(f"Run {int(r['run'])}: Observed fin temperature distribution")
        ax.set_ylabel("Temperature (��C)")
        ax.grid(True)

    axes[-1].set_xlabel("Position along fin (mm)")
    fig.tight_layout()
    fig.savefig( "observed_temperature_distribution.png", dpi=200)

    # ============================================================
    # PLOT 2: Theoretical temperature distribution (all tip models)
    # ============================================================
    fig, axes = plt.subplots(len(results), 1, figsize=(8, 12), sharex=True)
    if len(results) == 1:
        axes = [axes]

    x_dense_mm = np.linspace(0.0, l_fin_m, 300) * 1000.0

    for ax, (_, r) in zip(axes, results.iterrows()):
        obs_t = np.array([r["T1_C"], r["T2_C"], r["T3_C"], r["T4_C"], r["T5_C"]], dtype=float)
        t_amb = float(r["ambient_C"])
        tb = float(r["T1_C"])
        h_conv = float(r["h_W_m2K"])

        fin = fin_temperature_profiles(
            x_m=np.linspace(0.0, l_fin_m, 300),
            tb_c=tb,
            t_inf_c=t_amb,
            h=h_conv,
            k_fin=K_FIN_W_MK,
            d_fin_m=d_fin_m,
            l_fin_m=l_fin_m,
        )

        ax.plot(X_POS_MM, obs_t, "o", label="Observed")
        ax.plot(x_dense_mm, fin["t_ins"], label="Insulated tip")
        ax.plot(x_dense_mm, fin["t_ideal"], label="Ideal tip")
        ax.plot(x_dense_mm, fin["t_conv"], label="Convective tip")
        ax.set_title(f"Run {int(r['run'])}: Theoretical temperature distribution")
        ax.set_ylabel("Temperature (��C)")
        ax.grid(True)

    axes[-1].set_xlabel("Position along fin (mm)")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig("theoretical_temperature_distribution.png", dpi=200)

    # ============================================================
    # PLOT 3: Nu vs Re
    # ============================================================
    fig = plt.figure(figsize=(7, 5))
    plt.plot(results["Re"], results["Nu"], marker="o")
    plt.xlabel("Reynolds number (Re)")
    plt.ylabel("Nusselt number (Nu)")
    plt.title("Nu vs Re")
    plt.grid(True)
    fig.tight_layout()
    fig.savefig("nu_vs_re.png", dpi=200)

    # ============================================================
    # PLOT 4: Log-log plot of Nu vs Re
    # ============================================================
    fig = plt.figure(figsize=(7, 5))
    plt.loglog(results["Re"], results["Nu"], marker="o")
    plt.xlabel("Reynolds number (Re)")
    plt.ylabel("Nusselt number (Nu)")
    plt.title("Log-log plot of Nu vs Re")
    plt.grid(True, which="both")
    fig.tight_layout()
    fig.savefig( "nu_vs_re_loglog.png", dpi=200)

    # ============================================================
    # PLOT 5: Heat transfer rate comparison
    # ============================================================
    fig = plt.figure(figsize=(7, 5))
    plt.plot(results["run"], results["q_ins_W"], marker="o", label="Insulated tip")
    plt.plot(results["run"], results["q_ideal_W"], marker="o", label="Ideal tip")
    plt.plot(results["run"], results["q_conv_W"], marker="o", label="Convective tip")
    plt.xlabel("Run number")
    plt.ylabel("Heat transfer rate (W)")
    plt.title("Heat transfer rate through the fin")
    plt.grid(True)
    plt.legend()
    fig.tight_layout()
    fig.savefig("heat_transfer_rate_comparison.png", dpi=200)

    # ============================================================
    # PLOT 6: Effectiveness comparison
    # ============================================================
    fig = plt.figure(figsize=(7, 5))
    plt.plot(results["run"], results["effectiveness_ins"], marker="o", label="Insulated tip")
    plt.plot(results["run"], results["effectiveness_ideal"], marker="o", label="Ideal tip")
    plt.plot(results["run"], results["effectiveness_conv"], marker="o", label="Convective tip")
    plt.xlabel("Run number")
    plt.ylabel("Fin effectiveness")
    plt.title("Fin effectiveness comparison")
    plt.grid(True)
    plt.legend()
    fig.tight_layout()
    fig.savefig("effectiveness_comparison.png", dpi=200)



if __name__ == "__main__":
    main()
