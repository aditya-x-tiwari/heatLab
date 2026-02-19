import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================================
# CONSTANTS (From Sheet)
# ==========================================================

d = 0.025            # Diameter (m)
L = 0.176            # Test length (m)
k_water = 4180       # Cp of water (J/kgK)

# Cross-sectional area
A = np.pi * (d**2) / 4

# Positions along rod (9 points evenly spaced)
x = np.linspace(0, L, 9)

# ==========================================================
# READ CSV
# ==========================================================

df = pd.read_csv("input_data.csv")

V = df["Voltage_V"].values[0]
I = df["Current_A"].values[0]
m_dot = df["WaterFlow_kg_s"].values[0]

T_water_in = df["Water_In_C"].values[0]
T_water_out = df["Water_Out_C"].values[0]

# Rod temperatures
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
# CALCULATIONS
# ==========================================================

# Heat input
Q_elec = V * I

# Heat removed by water
Q_water = m_dot * k_water * (T_water_out - T_water_in)

# Average heat flow
Q_avg = (Q_elec + Q_water) / 2

# Temperature gradient (linear fit)
coeffs = np.polyfit(x, T, 1)
dT_dx = coeffs[0]

# Thermal conductivity
k = -Q_avg / (A * dT_dx)

# ==========================================================
# RADIAL LOSS (OPTIONAL INSIGHT)
# ==========================================================

Q_loss = Q_elec - Q_water

# ==========================================================
# SAVE RESULTS
# ==========================================================

results = pd.DataFrame({
    "Q_electrical_W": [Q_elec],
    "Q_water_W": [Q_water],
    "Q_average_W": [Q_avg],
    "Temperature_gradient_K_per_m": [dT_dx],
    "Thermal_conductivity_W_mK": [k],
    "Heat_loss_W": [Q_loss]
})

results.to_csv("rod_results.csv", index=False)

print("\n==============================")
print(f"Thermal Conductivity k = {k:.2f} W/mK")
print("==============================")

# ==========================================================
# PLOTS
# ==========================================================

# 1️⃣ Temperature Distribution
plt.figure()
plt.plot(x, T, marker='o')
plt.xlabel("Position along rod (m)")
plt.ylabel("Temperature (°C)")
plt.title("Temperature Distribution along Rod")
plt.grid(True)
plt.savefig("plot1.png")
plt.close()

# 2️⃣ Linear Fit Plot
plt.figure()
plt.plot(x, T, 'o', label="Experimental")
plt.plot(x, np.polyval(coeffs, x), linestyle='--', label="Linear Fit")
plt.xlabel("Position (m)")
plt.ylabel("Temperature (°C)")
plt.title("Linear Fit for Temperature Gradient")
plt.legend()
plt.grid(True)
plt.savefig("plot2.png")
plt.close()

# 3️⃣ Log Temperature Difference (advanced)
T_inf = T_water_in
theta = T - T_inf

plt.figure()
plt.semilogy(x, theta, marker='o')
plt.xlabel("Position (m)")
plt.ylabel("Temperature Excess (log scale)")
plt.title("Log Plot of Temperature Distribution")
plt.grid(True)
plt.savefig("plot3.png")
plt.close()
