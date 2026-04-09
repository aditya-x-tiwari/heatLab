import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================================
# CONSTANTS
# ==========================================================
d = 0.025                  # rod diameter (m)
L_total = 0.176            # total length (m)
cp_water = 4181            # J/kgK
K_ins = 0.34               # insulation conductivity (W/mK)

r_i = 0.04              # inner radius
r_o = 0.06              # outer radius 

A = np.pi * d**2 / 4
dx = L_total / 8           # spacing between 9 thermocouples

# ==========================================================
# READ DATA
# ==========================================================
df = pd.read_csv("input_data.csv")

m_dot = df["WaterFlow_kg_s"].values[0]
T_water_in = df["Water_In_C"].values[0]
T_water_out = df["Water_Out_C"].values[0]
T_AB_in = df["T_AB_in"].values[0]
T_AB_out = df["T_AB_out"].values[0]
T_BC_in = df["T_BC_in"].values[0]
T_BC_out = df["T_BC_out"].values[0]

T = np.array([
    df["T1"].values[0],
    df["T2"].values[0],
    df["T3"].values[0],
    df["T4"].values[0],
    df["T5"].values[0],
    df["T6"].values[0],
    df["T7"].values[0],
    df["T8"].values[0],
    df["T9"].values[0],
])

# ==========================================================
# SECTION AA
# ==========================================================
Q_AA = m_dot * cp_water * (T_water_out - T_water_in)
dTdx_AA = (T[8] - T[7]) / dx

k_AA = -Q_AA / (A * dTdx_AA)

T_mean_AA = (T[8] + T[7]) / 2

# ==========================================================
# RADIAL LOSS AB
# ==========================================================
L_AB = dx * 4  # adjust based on geometry

Q_radial_AB = (2*np.pi*K_ins*L_AB*(T_AB_out-T_AB_in)) / np.log(r_o/r_i)

Q_BB = Q_AA + Q_radial_AB

# ==========================================================
# SECTION BB
# ==========================================================
dTdx_BB = 0.5 * (
    (T[5] - T[4]) / dx +
    (T[4] - T[3]) / dx
)

k_BB = -Q_BB / (A * dTdx_BB)

T_mean_BB = (T[5] + T[4]) / 2

# ==========================================================
# RADIAL LOSS BC
# ==========================================================
L_BC = dx * 4

Q_radial_BC = (2*np.pi*K_ins*L_BC*(T_BC_out-T_BC_in)) / np.log(r_o/r_i)

Q_CC = Q_BB + Q_radial_BC

# ==========================================================
# SECTION CC
# ==========================================================
dTdx_CC = (T[1] - T[0]) / dx

k_CC = -Q_CC / (A * dTdx_CC)

T_mean_CC = (T[1] + T[0]) / 2

# ==========================================================
# PRINT RESULTS
# ==========================================================
print("\nThermal Conductivity Values:")
print(f"k_AA = {k_AA:.2f} W/mK")
print(f"k_BB = {k_BB:.2f} W/mK")
print(f"k_CC = {k_CC:.2f} W/mK")

# ==========================================================
# PLOTS
# ==========================================================
# Temperature vs Distance
x = np.linspace(0, L_total, 9)

# Plot for Temperature Distribution
plt.figure()
plt.scatter(x, T, marker='o', label="Experimental Data")  # Scatter plot (dots only)

# Linear fit for temperature vs distance (not necessary for this, but you can add if needed)
coeffs_temp = np.polyfit(x, T, 1)
fit_line_temp = np.poly1d(coeffs_temp)
plt.plot(x, fit_line_temp(x), 'r--', label=f'Best Fit: T = {coeffs_temp[0]:.2f}x + {coeffs_temp[1]:.2f}')

# Add equation to plot
equation_text_temp = f"T = {coeffs_temp[0]:.2f}x + {coeffs_temp[1]:.2f}"
plt.text(0.6, 85, equation_text_temp, transform=plt.gca().transAxes, fontsize=10, color='red')

plt.xlabel("Distance (m)")
plt.ylabel("Temperature (°C)")
plt.title("Temperature Distribution")
plt.grid(True)
plt.legend()
plt.savefig("T_vs_x.png")
plt.close()

# k vs Mean Temperature
k_values = [k_AA, k_BB, k_CC]
T_means = [T_mean_AA, T_mean_BB, T_mean_CC]

# Plot for k vs Mean Temperature
plt.figure()
plt.scatter(T_means, k_values, marker='o', label="Experimental Data")  # Scatter plot (dots only)

# Linear fit for k vs mean temperature
coeffs_k = np.polyfit(T_means, k_values, 1)
fit_line_k = np.poly1d(coeffs_k)
plt.plot(T_means, fit_line_k(T_means), 'r--', label=f'Best Fit: k = {coeffs_k[0]:.2f}T + {coeffs_k[1]:.2f}')

# Add equation to plot
equation_text_k = f"k = {coeffs_k[0]:.2f}T + {coeffs_k[1]:.2f}"
plt.text(0.6, 0.18, equation_text_k, transform=plt.gca().transAxes, fontsize=10, color='red')

plt.xlabel("Mean Temperature (°C)")
plt.ylabel("Thermal Conductivity (W/mK)")
plt.title("k vs Temperature")
plt.grid(True)
plt.legend()
plt.savefig("k_vs_T.png")
plt.close()
