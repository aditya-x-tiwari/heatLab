import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml

# ==========================================================
# CONSTANTS
# ==========================================================

Ri = 0.05   # Inner radius (m)
Ro = 0.10   # Outer radius (m)

# ==========================================================
# READ YAML
# ==========================================================

with open("input_data.yaml", "r") as file:
    yaml_data = yaml.safe_load(file)

exp = yaml_data["temperature_distribution_experiment"]

if not exp["enabled"]:
    raise ValueError("temperature_distribution_experiment is disabled in YAML")

data_list = exp["data"]

# ==========================================================
# CALCULATIONS
# ==========================================================

results = []

for data in data_list:

    V = data["Voltage_V"]
    I = data["Current_A"]

    T = np.array(data["T_values"])

    # Inner temps → first 4
    Ti = np.mean(T[:4])

    # Outer temps → remaining
    To = np.mean(T[4:])

    # Heat input
    Q = V * I

    # Thermal conductivity
    k = (Q * (Ro - Ri)) / (4 * np.pi * Ri * Ro * (Ti - To))

    results.append([V, I, Ti, To, Q, k])

# ==========================================================
# SAVE RESULTS
# ==========================================================

results_df = pd.DataFrame(results, columns=[
    "Voltage_V", "Current_A",
    "T_inner_avg_C", "T_outer_avg_C",
    "Heat_Input_W", "Thermal_Conductivity_W_mK"
])

results_df.to_csv("spherical_results.csv", index=False)

print("\n==== RESULTS ====\n")
print(results_df)

# ==========================================================
# PLOTS
# ==========================================================

# 1️⃣ k vs Heat Input
plt.figure()
plt.plot(results_df["Heat_Input_W"], results_df["Thermal_Conductivity_W_mK"], marker='o')
plt.xlabel("Heat Input (W)")
plt.ylabel("Thermal Conductivity (W/mK)")
plt.title("k vs Heat Input")
plt.grid(True)
plt.show()

# 2️⃣ k vs Temperature Difference
delta_T = results_df["T_inner_avg_C"] - results_df["T_outer_avg_C"]

plt.figure()
plt.plot(delta_T, results_df["Thermal_Conductivity_W_mK"], marker='o')
plt.xlabel("Temperature Difference (Ti - To)")
plt.ylabel("Thermal Conductivity (W/mK)")
plt.title("k vs Temperature Difference")
plt.grid(True)
plt.show()

# 3️⃣ Voltage vs Heat Input
plt.figure()
plt.plot(results_df["Voltage_V"], results_df["Heat_Input_W"], marker='o')
plt.xlabel("Voltage (V)")
plt.ylabel("Heat Input (W)")
plt.title("Voltage vs Heat Input")
plt.grid(True)
plt.show()
