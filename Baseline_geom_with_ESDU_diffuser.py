from modules import geometry
from modules import esdu_75026 as esdu
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 20})

# sigma = np.array([
#     [1,     1.3],
#     [1.1,   1.4],
#     [1.2,   1.5]
# ])
sigma = np.array([
    [1,     1.2],
    [1.1,   1.3],
])
length_factor = 0.2


nrows, ncols = sigma.shape
width = 6 * ncols
height = 4 * nrows

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(width, height), sharex='all', sharey='all')

for i in range(nrows):
    for j in range(ncols):
        sigma_ij = sigma[i, j]

        Cpr2, L_h1 = esdu.optimum_Cpr2_symmetrical(sigma_ij)
        L_h1 *= length_factor
        angle = esdu.calc_phi(sigma_ij, L_h1)

        geom = geometry.NonDimensionalGeometry(
            Lh_intake=0.4,   
            Lh_rotor=0.4,     
            Lh_diffuser=L_h1,    
            rr_hub_tip=0.3, 
            th_cowl=0.88,
            tc_rotor=0.1,
            tc_stator=0.1,
            phi_i=angle,
            phi_o=angle,
            cL_rotor=0.6,
            cL_stator=0.9,
            delta_oh=0.02,
            N_rotor=6,
            N_stator=7
        )

        xr_hub = geom.calc_hub_line()
        xr_cas = geom.calc_cas_line()
        xr_cowl = geom.calc_cowl_line()
        curves = [xr_hub, xr_cas, xr_cowl]

        for c in curves:
            axs[i, j].plot(c[0, :], c[1, :], 'r-')

        L_D = geom.Lh_total / (2 * geom.Rh_cas)
        axs[i, j].set_title(fr"$ \sigma = {sigma_ij} $, $ L/D = {L_D:.2f}$")
        axs[i, j].set_aspect('equal')
        axs[i, j].grid('major')

        if i == nrows-1:
            axs[i, j].set_xlabel(r"$ \dfrac{x}{h_r} $")
        if j == 0:
            axs[i, j].set_ylabel(r"$ \dfrac{r}{h_r} $", rotation=0)

fig.set_constrained_layout(True)

if length_factor == 1:
    fig.savefig("Figures/Baseline_geom_with_ESDU_diffuser.pdf")
else:
    fig.savefig(f"Figures/Baseline_geom_with_ESDU_diffuser_{length_factor}fac.pdf")

plt.show()

