import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def process_pin_fin(csv_file):
    df = pd.read_csv(csv_file)

    # ---------------- CONSTANTS ----------------
    g = 9.81
    rho_air = 1.2          # kg/m3
    mu_air = 1.8e-5        # Pa.s
    k_air = 0.026          # W/mK
    k_fin = 200            # W/mK (Al approx)
    Cd = 0.62

    D_fin = 0.0127         # m
    L = 0.15               # m
    D_orifice = 0.01       # m
    A_orifice = np.pi * D_orifice**2 / 4
    A_duct = 0.1 * 0.15    # m2

    # ---------------- AIR VELOCITY ----------------
    # Manometer head in meters
    df['h_man'] = df['manometer_diff_cm'] / 100

    # Orifice theoretical velocity
    df['V'] = Cd * np.sqrt(2 * g * df['h_man'])

    # ---------------- REYNOLDS ----------------
    df['Re'] = (rho_air * df['V'] * D_fin) / mu_air

    # ---------------- NUSSELT ----------------
    df['Nu'] = np.where(df['Re'] < 4000,
                        0.615 * df['Re']**0.466,
                        0.174 * df['Re']**0.618)

    # ---------------- HEAT TRANSFER COEFF ----------------
    df['h'] = df['Nu'] * k_air / D_fin

    # ---------------- FIN PARAMETERS ----------------
    P = np.pi * D_fin
    A = np.pi * D_fin**2 / 4

    df['m'] = np.sqrt(df['h'] * P / (k_fin * A))

    # Base temperature assumed T1
    df['Tb'] = df['T1']
    df['Tinf'] = df['T_ambient']

    df['Q_fin'] = np.sqrt(df['h'] * P * k_fin * A) * \
                  (df['Tb'] - df['Tinf']) * \
                  np.tanh(df['m'] * L)

    df['Efficiency'] = np.tanh(df['m'] * L) / (df['m'] * L)

    print(df[['Re','h','Q_fin','Efficiency']])

    # -------- PLOTS --------

    # Temperature distribution
    x_positions = np.linspace(0, L, 5)

    for i, row in df.iterrows():
        temps = [row['T1'], row['T2'], row['T3'],
                 row['T4'], row['T5']]
        plt.plot(x_positions, temps)

    plt.xlabel("Fin Length (m)")
    plt.ylabel("Temperature (°C)")
    plt.title("Temperature Distribution Along Fin")
    plt.grid(True)
    plt.savefig("Temperature_vs_Length.png")
    plt.close()

    # h vs Re
    plt.figure()
    plt.plot(df['Re'], df['h'], 'o-')
    plt.xlabel("Reynolds Number")
    plt.ylabel("Heat Transfer Coefficient (h)")
    plt.title("h vs Re")
    plt.grid(True)
    plt.savefig("h_vs_Re.png")
    plt.close()

    # Q vs Velocity
    plt.figure()
    plt.plot(df['V'], df['Q_fin'], 'o-')
    plt.xlabel("Air Velocity (m/s)")
    plt.ylabel("Heat Transfer Rate (W)")
    plt.title("Heat Transfer vs Velocity")
    plt.grid(True)
    plt.savefig("Q_vs_Velocity.png")
    plt.close()

    df.to_csv("pin_fin_results.csv", index=False)

    print("Files saved:")
    print(" - Temperature_vs_Length.png")
    print(" - h_vs_Re.png")
    print(" - Q_vs_Velocity.png")
    print(" - pin_fin_results.csv")


if __name__ == "__main__":
    process_pin_fin("pin_fin_input.csv")
