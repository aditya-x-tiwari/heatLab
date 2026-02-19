import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml

# ==========================================================
# CONSTANTS (From Sheet)
# ==========================================================

d = 0.025            # Diameter (m)
L = 0.176            # Test length (m)
k_water = 4180       # Cp of water (J/kgK)

A = np.pi * (d**2) / 4
x = np.linspace(0, L, 9)

# ==========================================================
# READ YAML
# ==========================================================

with open("input_data.yaml", "r") as file:
    yaml_data = yaml.safe_load(file)

exp = yaml_data["water_flow_experiment"]

if not exp["enabled"]:
    raise ValueError("water_flow_experiment is disabled in YAML")

data = exp["data"]

# Extract values
V = data["Voltage_V"]
I = data["Current_A"]
m_dot = data["WaterFlow_kg_s"]

T_water_in = data["Water_In_C"]
T_water_out = data["Water_Out_C"]

# Temperatures array
T = np.array(data["T_values"])

# ==========================================================
# CALCULATIONS
# ==========================================================

Q_elec = V * I
Q_water = m_dot * k_water * (T_water_out - T_water_in)
Q_avg = (Q_elec + Q_water) / 2

coeffs = np.polyfit(x, T, 1)
dT_dx = coeffs[0]

k = -Q_avg / (A * dT_dx)

# ==========================================================
# RADIAL LOSS
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

plt.figure()
plt.plot(x, T, marker='o')
plt.xlabel("Position along rod (m)")
plt.ylabel("Temperature (°C)")
plt.title("Temperature Distribution along Rod")
plt.grid(True)
plt.show()

plt.figure()
plt.plot(x, T, 'o', label="Experimental")
plt.plot(x, np.polyval(coeffs, x), linestyle='--', label="Linear Fit")
plt.xlabel("Position (m)")
plt.ylabel("Temperature (°C)")
plt.title("Linear Fit for Temperature Gradient")
plt.legend()
plt.grid(True)
plt.show()

T_inf = T_water_in
theta = T - T_inf

plt.figure()
plt.semilogy(x, theta, marker='o')
plt.xlabel("Position (m)")
plt.ylabel("Temperature Excess (log scale)")
plt.title("Log Plot of Temperature Distribution")
plt.grid(True)
plt.show()

