import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI


# ==========================================================
# READ CSV
# ==========================================================

df = pd.read_csv("input_data.csv")

H_mm = df["H_mm"].values
V = df["Voltage_V"].values
Ts = df["Te_C"].values + 273.15
Ta = df["Td_C"].values + 273.15


# ==========================================================
# CONSTANTS (From Lab Sheet)
# ==========================================================

duct = 65*150e-6        # Duct Cross-Section (m^2)
d = 0.0158              # Diameter (m)
l = 0.05                # Length (m)
A = 2.482e-3            # Heated area (m^2)
R = 70                  # Resistance (ohm)
C = 74.294              # Velocity constant from sheet

Pa = 775*13600*9.81/1000 # Atmospheric Pressure (Pa) 
Tf = (Ts + Ta) / 2      # Film temperature

rho_s = PropsSI("D", "T", Ts, "P", Pa, "Air")
mu_s  = PropsSI("VISCOSITY", "T", Ts, "P", Pa, "Air")
nu_s = mu_s / rho_s           # Kinematic viscosity

rho_a = PropsSI("D", "T", Ta, "P", Pa, "Air")
mu_a  = PropsSI("VISCOSITY", "T", Ta, "P", Pa, "Air")
nu_a = mu_a / rho_a           # Kinematic viscosity

k_a   = PropsSI("CONDUCTIVITY", "T", Ta, "P", Pa, "Air")   # Thermal conductivity of air (W/m.K)
Pr_a  = PropsSI("PRANDTL", "T", Ta, "P", Pa, "Air")        # Prandtl's Number of air 



# ==========================================================
# CALCULATIONS
# ==========================================================

# Convert mm to meters
H = H_mm / 1000

# Air velocity
U = C * np.sqrt(H*Ta/Pa)

# Electrical heat input
Q = (V**2) / R

# Heat flux
q = Q / A

# Reynolds number at film temperature
rho_f = PropsSI("D", "T", Tf, "P", Pa, "Air")
mu_f  = PropsSI("VISCOSITY", "T", Tf, "P", Pa, "Air")
nu_f = mu_f / rho_f          # Kinematic viscosity
Re_f = (U * d) / nu_f

# Reynolds number
Re_a = (U * d) / nu_a

# Heat transfer coefficient
h = q / (Ts - Ta)

# Nusselt number
Nu = (h * d) / k_a

# ==========================================================
# CORRELATIONS (From Sheet)
# ==========================================================


flow_regime = []

for r in Re_f:
    if r < 4000:
        Nu_1 = 0.615 * (r ** 0.466)
    elif 4000 <= r < 40000:
        Nu_1 = 0.174 * (r ** 0.618)
    else:
        Nu_1 = 0.0239 * (r ** 0.805)

for r in Re_a:
        Nu_2 = (0.4 * (r ** 0.5) + 0.06 * (r ** 0.667) ) * (Pr_a ** 0.4) * ((mu_a/mu_s) ** 0.25)



# ==========================================================
# SAVE RESULTS
# ==========================================================


df["Heat_Input_W"] = Q
df["Heat_Flux_W_m2"] = q
df["Heat_Transfer Coefficient_W/m2 k"] = h
df["Velocity_m_s"] = U
df["Re_a"] = Re_a
df["Re_f"] = Re_f
df["Nu"] = Nu
df["Nu_corr_1"] = Nu_1
df["Nu_corr_2"] = Nu_2

df.to_csv("output_results.csv", index=False)

print("\nResults saved to output_results.csv")

# ==========================================================
# PLOTS
# ==========================================================



#  Log-Log Plot of Nu vs Re_a 
plt.figure()
plt.loglog(Re_a, Nu, marker='o', label="Experimental")
plt.xlabel("Re_a")
plt.ylabel("Nu")
plt.title("Log-Log: Nu vs Re")
plt.legend()
plt.grid(True, which="both")
plt.savefig("plot1.png")
plt.close()


#  Log-Log Plot of Nu_1 vs Re_f
plt.figure()
plt.loglog(Re_f, Nu_1, marker='o', label="Correlation 1")
plt.xlabel("Re_f")
plt.ylabel("Nu_1")
plt.title("Log-Log: Nu_1 vs Re_f")
plt.legend()
plt.grid(True, which="both")
plt.savefig("plot2.png")
plt.close()


#  Log-Log Plot of Nu_2 vs Re_a 
plt.figure()
plt.loglog(Re_a, Nu_2, marker='o', label="Correlation 2")
plt.xlabel("Re_a")
plt.ylabel("Nu_2")
plt.title("Log-Log: Nu_2 vs Re_a")
plt.legend()
plt.grid(True, which="both")
plt.savefig("plot3.png")
plt.close()

