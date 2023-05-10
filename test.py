from geometry import *


geom = Geometry.example()

xr_hub = calc_hub_line(geom)
xr_cas = calc_cas_line(geom)
xr_cowl = calc_cowl_line(geom)

plt.plot(xr_hub[0, :], xr_hub[1, :])
plt.plot(xr_cas[0, :], xr_cas[1, :])
plt.plot(xr_cowl[0, :], xr_cowl[1, :])
plt.axis("equal")
plt.show()