# Two-Compartment PK Model with Non-linear ODE Solver

## 📌 Project Overview
This project implements a mechanistic Pharmacokinetic (PK) model representing a **Two-Compartment system with non-linear (saturable) elimination kinetics**. 

Instead of assuming standard first-order elimination, this architecture incorporates Michaelis-Menten kinetics to model metabolic pathways that undergo enzyme saturation at higher concentrations (akin to saturation phenomena modeled in chemical-biological systems).

## 📊 Mathematical Architecture
The temporal distribution and elimination profiles of the drug are governed by the following system of non-linear Ordinary Differential Equations (ODEs):

$$ \frac{dAc}{dt} = - Q \cdot \left(\frac{Ac}{Vc} - \frac{Ap}{Vp}\right) - \frac{V_{max} \cdot \frac{Ac}{Vc}}{K_m + \frac{Ac}{Vc}} $$

$$ \frac{dAp}{dt} = Q \cdot \left(\frac{Ac}{Vc} - \frac{Ap}{Vp}\right) $$

Where:
* $Ac, Ap$: Amount of drug in central and peripheral compartments (mg).
* $Vc, Vp$: Apparent volumes of distribution (L).
* $Q$: Intercompartmental clearance (L/h).
* $V_{max}$: Maximum metabolic elimination capacity (mg/h).
* $K_m$: Michaelis constant (mg/L).

## 🛠️ Implementation Details
* **Solver Selection:** Implemented using `scipy.integrate.solve_ivp` with the **Radau** method, ensuring mathematical convergence and absolute stability even under stiff conditions caused by rapid saturation transitions.
* **Outputs:** The script generates time-concentration profiles for both compartments and shows a high-resolution simulation plot.
  
  <img width="781" height="468" alt="download" src="https://github.com/user-attachments/assets/780d5371-1e1d-49d4-a86b-31ce745909f6" />


## 🚀 How to Run
```bash
pip install numpy scipy matplotlib
python pk_solver.py
