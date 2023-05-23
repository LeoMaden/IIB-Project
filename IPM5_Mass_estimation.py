import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from matplotlib.ticker import MultipleLocator

from modules import geometry, propulsor_mass

plt.rcParams.update({"font.size": 20})

# Load IPM5 geometry
curves = scipy.io.loadmat("Fan Curves/IPM5_Curves.mat")["c"]
splines = scipy.io.loadmat("Fan Curves/IPM5_splines.mat", squeeze_me=True)["B"]

# Calculations
Rc = splines[0]["xr_cas"].item()[0, 1]
Rh = splines[0]["xr_hub"].item()[0, 1]
hr = Rc - Rh
D = 2 * Rc
A3 = np.pi * (Rc**2 - Rh**2)
IPM5_mass = 2      # Actual IPM5 mass in kg

# Plots
fig1, ax1 = plt.subplots(figsize=(10, 8))
fig2, ax2 = plt.subplots(figsize=(10, 8))

# Plot IPM5 geometry
c_names = ["xr_cowl", "xr_int", "xr_nose", "xr_exit", "xr_tail"]

for name in c_names:
    xr = curves[name].item() / hr
    ax1.plot(xr[:, 0], xr[:, 1], 'k')

for i in [0, 1]:
    for name in ["xr_hub", "xr_cas"]:
        xr = splines[i][name].item() / hr
        
        ax1.plot(xr[:, 0], xr[:, 1], 'k')

ax1.set_aspect("equal")

black_line = mlines.Line2D([], [], color='k', linestyle='-')
red_line = mlines.Line2D([], [], color='r', linestyle='-')
ax1.legend(handles=[black_line, red_line], labels=["Actual", "Model"])

ax1.set_xlabel(r"$ \dfrac{x}{h_r} $")
ax1.set_ylabel(r"$ \dfrac{r}{h_r} $", rotation=0)

fig1.tight_layout()

# Generate model geometry
geom = geometry.NonDimensionalGeometry(
            Lh_intake=0.7,   
            Lh_rotor=1.1,     
            Lh_diffuser=1.45,    
            rr_hub_tip=0.35, 
            th_cowl=0.42,
            tc_rotor=0.1,
            tc_stator=0.1,
            phi_i=20,
            phi_o=-3,
            cL_rotor=0.2,
            cL_stator=0.2,
            delta_oh=0.02,
            N_rotor=9,
            N_stator=10
        )

xr_hub = geom.calc_hub_line()
xr_cas = geom.calc_cas_line()
xr_cowl = geom.calc_cowl_line()

curves = [xr_hub, xr_cas, xr_cowl]

# Plot model geometry
for c in curves:
    ax1.plot(c[0, :], c[1, :], 'r')

# Mass parameters
rho_air = 1.2
rho_nylon = 1150
rho_cu = 8960
rho_al = 2700
rho_resin = 1200

mass_params = propulsor_mass.NonDimensionalMassParams(
    cowl_density        =   rho_nylon / rho_air,
    hub_density         =   rho_cu / rho_air,
    rotor_density       =   rho_al / rho_air,
    stator_density      =   rho_resin / rho_air,
    cowl_solidity       =   0.2,
    hub_solidity        =   0.5
)

# Calculate nondimensional mass
mass_model = propulsor_mass.NonDimensionalMassModel(mass_params, geom)

rotor_mass_nd = mass_model.calc_M_rotor()
stator_mass_nd = mass_model.calc_M_stator()
hub_mass_nd = mass_model.calc_M_hub(xr_hub)
cowl_mass_nd = mass_model.calc_M_cowl(xr_cas, xr_cowl)

# Dimensionalise mass
norm = rho_air * A3 * D

rotor_mass = norm * rotor_mass_nd
stator_mass = norm * stator_mass_nd
hub_mass = norm * hub_mass_nd
cowl_mass = norm * cowl_mass_nd

total_mass = rotor_mass + stator_mass + hub_mass + cowl_mass

print(f"Estimated mass = {total_mass}")

# Plot mass breakdown on bar chart
ax2.bar("Model", rotor_mass, label='Rotor')
ax2.bar("Model", stator_mass, bottom=rotor_mass, label='Stator')
ax2.bar("Model", hub_mass, bottom=rotor_mass+stator_mass, label='Hub')
ax2.bar("Model", cowl_mass, bottom=rotor_mass+stator_mass+hub_mass, label='Cowl')

ax2.bar("Actual", IPM5_mass, color='grey')

ax2.legend(bbox_to_anchor=(0, 1.001, 1, 0.1), loc='lower left', ncol=4, mode="expand")
ax2.set_ylabel("Mass (kg)")
ax2.yaxis.set_major_locator(MultipleLocator(0.5))
ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
fig2.tight_layout()

fig1.savefig("Figures/IPM5_Model_vs_actual_geom.pdf")
fig2.savefig("Figures/IPM5_Model_vs_actual_mass.pdf")

plt.show()