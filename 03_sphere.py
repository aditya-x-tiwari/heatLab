import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================================
# CONSTANTS
# ==========================================================

Ri = 0.05   # Inner radius (m)
Ro = 0.10   # Outer radius (m)

# ==========================================================
# READ CSV
# ==========================================================

df = pd.read_csv("input_data.csv")

# ==========================================================
# CALCULATIONS
# ==========================================================

results = []

for index, row in df.iterrows():

    V = row["Voltage_V"]
    I = row["Current_A"]

    # Inner temps (T1–T4)
    Ti = np.mean([row["T1"], row["T2"], row["T3"], row["T4"]])

    # Outer temps (T5–T10)
    To = np.mean([
        row["T5"], row["T6"], row["T7"],
        row["T8"], row["T9"], row["T10"]
    ])

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
plt.savefig("plot1.png")
plt.close()

# 2️⃣ k vs Inner Surface Temperature 
plt.figure()
plt.plot(results_df["T_inner_avg_C"] , results_df["Thermal_Conductivity_W_mK"], marker='o')
plt.xlabel("Inner Surface Temperature (Ti)")
plt.ylabel("Thermal Conductivity (W/mK)")
plt.title("k vs Inner Surface Temperature ")
plt.grid(True)
plt.savefig("plot2.png")
plt.close()

