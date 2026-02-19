import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================================
# CONSTANTS (From Lab Sheet)
# ==========================================================

duct = 65*150e-6        # Duct Cross-Section (m^2)
d = 0.0158              # Diameter (m)
l = 0.05                # Length (m)
A = 2.482e-6            # Heated area (m^2)
R = 70                  # Resistance (ohm)
C = 74.294              # Velocity constant from sheet
k = 0.026               # Thermal conductivity of air (W/m.K)
rho = 1.165             # Density of air (kg/m3 at ~30C)
mu = 1.85e-5            # Dynamic viscosity (Pa.s at ~30C)
Pa = 775*13600*9.81/1000 # Atmospheric Pressure (Pa) 

nu = mu / rho           # Kinematic viscosity

# ==========================================================
# READ CSV
# ==========================================================

df = pd.read_csv("input_data.csv")

H_mm = df["H_mm"].values
V = df["Voltage_V"].values
Te = df["Te_C"].values
Td = df["Td_C"].values

# ==========================================================
# CALCULATIONS
# ==========================================================

# Convert mm to meters
H = H_mm / 1000

# Air velocity
U = C * np.sqrt(H*Te/Pa)

# Electrical heat input
Q = (V**2) / R

# Heat flux
q = Q / A

# Film temperature
Tf = (Te + Td) / 2

# Reynolds number
Re = (U * d) / nu

# Heat transfer coefficient
h = q / (Te - Td)

# Nusselt number
Nu = (h * d) / k

# ==========================================================
# FLOW REGIME
# ==========================================================

flow_regime = []

for r in Re:
    if r < 4000:
        flow_regime.append("Laminar")
    elif 4000 <= r < 40000:
        flow_regime.append("Transition")
    else:
        flow_regime.append("Turbulent")

# ==========================================================
# CORRELATIONS (From Sheet)
# ==========================================================

Nu_corr_1 = 0.615 * (Re ** 0.466)
Nu_corr_2 = 0.174 * (Re ** 0.618)
Nu_corr_3 = 0.0239 * (Re ** 0.805)

# ==========================================================
# POWER LAW FIT (Experimental)
# ==========================================================

log_Re = np.log(Re)
log_Nu = np.log(Nu)

n, log_C = np.polyfit(log_Re, log_Nu, 1)
C_exp = np.exp(log_C)

Nu_fit = C_exp * (Re ** n)

# ==========================================================
# SAVE RESULTS
# ==========================================================

df["Velocity_m_s"] = U
df["Heat_Input_W"] = Q
df["Heat_Flux_W_m2"] = q
df["Film_Temp_C"] = Tf
df["Re"] = Re
df["Flow_Regime"] = flow_regime
df["h_W_m2K"] = h
df["Nu"] = Nu
df["Nu_corr_1"] = Nu_corr_1
df["Nu_corr_2"] = Nu_corr_2
df["Nu_corr_3"] = Nu_corr_3
df["Nu_fit"] = Nu_fit

df.to_csv("output_results.csv", index=False)

print("\n==============================")
print("Experimental Power Law Fit:")
print(f"Nu = {C_exp:.4f} Re^{n:.4f}")
print("==============================")
print("\nResults saved to output_results.csv")

# ==========================================================
# PLOTS
# ==========================================================

# 1️⃣ Nu vs Re
plt.figure()
plt.plot(Re, Nu, marker='o', label="Experimental")
plt.plot(Re, Nu_corr_1, linestyle='--', label="Corr 1")
plt.plot(Re, Nu_corr_2, linestyle='--', label="Corr 2")
plt.plot(Re, Nu_corr_3, linestyle='--', label="Corr 3")
plt.xlabel("Re")
plt.ylabel("Nu")
plt.title("Nu vs Re")
plt.legend()
plt.grid(True)
plt.savefig("plot1.png")
plt.close()

# 2️⃣ Log-Log Plot
plt.figure()
plt.loglog(Re, Nu, marker='o', label="Experimental")
plt.loglog(Re, Nu_fit, linestyle='--', label="Power Law Fit")
plt.xlabel("Re")
plt.ylabel("Nu")
plt.title("Log-Log: Nu vs Re")
plt.legend()
plt.grid(True, which="both")
plt.savefig("plot2.png")
plt.close()

# 3️⃣ h vs Re
plt.figure()
plt.plot(Re, h, marker='o')
plt.xlabel("Re")
plt.ylabel("h (W/m²K)")
plt.title("h vs Re")
plt.grid(True)
plt.savefig("plot3.png")
plt.close()

# 4️⃣ Velocity vs Re
plt.figure()
plt.plot(Re, U, marker='o')
plt.xlabel("Re")
plt.ylabel("Velocity (m/s)")
plt.title("Velocity vs Re")
plt.grid(True)
plt.savefig("plot4.png")
plt.close()
