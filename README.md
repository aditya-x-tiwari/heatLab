# heatLab 🔥

Computational toolkit for Heat Transfer Laboratory experiments.

This repository contains Python-based automation scripts developed for analyzing heat transfer experiments, including forced convection over pin fins.

Developed as part of undergraduate Mechanical Engineering coursework.

---

## 📌 Experiments Covered

### 1️⃣ Heat Transfer from a Pin Fin (Forced Convection)

Includes computation of:

- Air velocity from orifice manometer
- Reynolds number
- Nusselt number
- Convective heat transfer coefficient (h)
- Fin parameter (m)
- Heat transfer rate (Q_fin)
- Fin efficiency (η)

---

## 📊 Generated Outputs

- Temperature distribution along fin (T vs x)
- h vs Reynolds number
- Heat transfer rate vs velocity
- Computed performance metrics exported to CSV

---

## 🔬 Implemented Theory

Reynolds Number:
Re = ρVD / μ

Nusselt Correlation (Forced Convection):
Nu = 0.615 Re^0.466  (40 < Re < 4000)

Heat Transfer Coefficient:
h = Nu·k / D

Fin Efficiency:
η = tanh(mL) / (mL)

Heat Transfer Rate:
Q = √(hPkA) (Tb − T∞) tanh(mL)

---

## ⚙️ Tech Stack

- Python 3.x
- NumPy
- Pandas
- Matplotlib

---

## 📎 How to Use

1. Input experimental readings into CSV.
2. Run the script:
   ```bash
   python {experiment_name}.py
