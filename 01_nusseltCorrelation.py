import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# CONSTANTS (From Experiment Sheet)
# -----------------------------
d = 0.0158              # Diameter (m)
A = 2.482e-2            # Heated area (m^2)
R = 70                  # Resistance (ohm)
C = 4.294               # Velocity constant from sheet
k = 0.026               # Thermal conductivity of air (W/m.K)
rho = 1.165             # Density of air (kg/m3)
mu = 1.85e-5            # Dynamic viscosity (Pa.s)

# -----------------------------
# READ CSV FILE
# -----------------------------
df = pd.read_csv("input_data.csv")

# Extract columns
H_mm = df["H_mm"].values
V = df["Voltage_V"].values
Te = df["Te_C"].values
Td = df["Td_C"].values

# -----------------------------
# CALCULATIONS
# -----------------------------

# Convert H to meters
H = H_mm / 1000

# Air velocity
U = C * np.sqrt(H)

# Heat input
Q = (V**2) / R

# Heat flux
q = Q / A

# Film temperature
Tf = (Te + Td) / 2

# Kinematic viscosity
nu = mu / rho

# Reynolds number
Re = (U * d) / nu

# Heat transfer coefficient
h = q / (Te - Td)

# Nusselt number
Nu = (h * d) / k

# -----------------------------
# CORRELATION VALUES (From Sheet)
# -----------------------------
Nu_corr_1 = 0.615 * (Re**0.466)
Nu_corr_2 = 0.174 * (Re**0.618)
Nu_corr_3 = 0.0239 * (Re**0.805)

# -----------------------------
# SAVE RESULTS
# -----------------------------
df["Velocity_m_s"] = U
df["Heat_Input_W"] = Q
df["Re"] = Re
df["h_W_m2K"] = h
df["Nu"] = Nu
df["Nu_corr_1"] = Nu_corr_1
df["Nu_corr_2"] = Nu_corr_2
df["Nu_corr_3"] = Nu_corr_3

df.to_csv("output_results.csv", index=False)

print("\nAll values computed and saved to output_results.csv")

# -----------------------------
# PLOTS
# -----------------------------

# 1. Nu vs Re
plt.figure()
plt.plot(Re, Nu, marker='o')
plt.xlabel("Re")
plt.ylabel("Nu")
plt.title("Nu vs Re")
plt.show()

# 2. Log-Log Plot
plt.figure()
plt.loglog(Re, Nu, marker='o')
plt.xlabel("Re")
plt.ylabel("Nu")
plt.title("Log-Log: Nu vs Re")
plt.show()

# 3. h vs Re
plt.figure()
plt.plot(Re, h, marker='o')
plt.xlabel("Re")
plt.ylabel("h (W/m2K)")
plt.title("h vs Re")
plt.show()

# 4. Velocity vs Re
plt.figure()
plt.plot(Re, U, marker='o')
plt.xlabel("Re")
plt.ylabel("Velocity (m/s)")
plt.title("Velocity vs Re")
plt.show()
