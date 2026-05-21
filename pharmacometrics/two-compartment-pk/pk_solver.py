import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def twocomp_pk_system(t, y, Vmax, Km, Q, Vc, Vp):
    """
    Defines the differential equations for a 2-compartment PK model
    with non-linear (Michaelis-Menten) elimination.
    """
    Ac, Ap = y  # Amounts in Central and Peripheral compartments (mg)
    
    # Concentrations
    Cc = Ac / Vc  # Central concentration (mg/L)
    Cp = Ap / Vp  # Peripheral concentration (mg/L)
    
    # Distribution fluxes between compartments
    flux_distribution = Q * (Cc - Cp)
    
    # Non-linear elimination flux (Saturable metabolism)
    flux_elimination = (Vmax * Cc) / (Km + Cc)
    
    # ODEs
    dAc_dt = - flux_distribution - flux_elimination
    dAp_dt = flux_distribution
    
    return [dAc_dt, dAp_dt]

def run_simulation():
    # --- Simulated Parametric Space (Anonymized) ---
    Vc = 5.0      # Volume of central compartment (L)
    Vp = 10.0     # Volume of peripheral compartment (L)
    Q = 0.8       # Intercompartmental clearance (L/h)
    Vmax = 12.0   # Max elimination rate (mg/h)
    Km = 2.5      # Concentration at half-Vmax (mg/L)
    
    # Initial conditions: IV Bolus of 100 mg injected into Central Compartment
    Initial_Dose = 100.0 # mg
    y0 = [Initial_Dose, 0.0]
    
    # Time vector (0 to 24 hours)
    t_span = (0, 24)
    t_eval = np.linspace(0, 24, 200)
    
    # --- Numerical Solution ---
    sol = solve_ivp(
        twocomp_pk_system, 
        t_span, 
        y0, 
        args=(Vmax, Km, Q, Vc, Vp), 
        t_eval=t_eval,
        method='Radau' # Stiff-tolerant solver to ensure mathematical robustness
    )
    
    # --- Fit & Visualization Diagnostics ---
    time = sol.t
    central_conc = sol.y[0] / Vc
    peripheral_conc = sol.y[1] / Vp
    
    plt.figure(figsize=(9, 5))
    plt.plot(time, central_conc, label='Central Compartment (Plasma)', color='navy', lw=2)
    plt.plot(time, peripheral_conc, label='Peripheral Compartment (Tissue)', color='crimson', lw=2, linestyle='--')
    plt.title('Two-Compartment PK Profile with Saturable (Non-Linear) Elimination', fontsize=12, fontweight='bold')
    plt.xlabel('Time (Hours)', fontsize=10)
    plt.ylabel('Concentration (mg/L)', fontsize=10)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(frameon=True)
    
    # Save plot as image for the documentation
plt.show()

if __name__ == "__main__":
    run_simulation()


