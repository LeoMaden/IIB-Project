from geometry import *
import esdu_75026 as esdu







# --- ESDU Test ---
# A2_A1 = 3
# Cpr2, L_h1 = esdu.optimum_Cpr2_symmetrical(A2_A1)
# print(f"Area ratio = {A2_A1}")
# print(f"Optimum Cpr = {Cpr2}")
# print(f"When diffuser is length L/h1 = {L_h1}")



# --- Geom test ---
# geom = NonDimensionalGeometry.example()
# xr_hub = geom.calc_hub_line()
# xr_cas = geom.calc_cas_line()
# xr_cowl = geom.calc_cowl_line()

# print("Hub volume = {:.3g}".format(geom.calc_Vh3_hub(xr_hub)))
# print("Cowl volume = {:.3g}".format(geom.calc_Vh3_cowl(xr_cas, xr_cowl)))

# plt.plot(xr_hub[0, :], xr_hub[1, :], marker="")
# plt.plot(xr_cas[0, :], xr_cas[1, :], marker="")
# plt.plot(xr_cowl[0, :], xr_cowl[1, :], marker="")
# #plt.axis("off")
# plt.axis("equal")
# plt.show()