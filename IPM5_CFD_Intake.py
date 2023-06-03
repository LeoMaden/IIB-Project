import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import matplotlib.colors as mcolors

plt.rcParams.update({"font.size": 20})

rpm = 22000

# Air viscosity
mu = 1.81e-5

# Design params
D = 0.16
hub_tip = 0.35
phi_des = 0.7

Rc = D / 2
Rh = Rc * hub_tip
Rmean = 0.5 * (Rc + Rh)


# Load data
filename = f"CFD Data/IPM5_intake_{rpm}rpm.mat"
data = scipy.io.loadmat(filename, squeeze_me=True, chars_as_strings=True)

# alpha, Vf, i, k
x = data['x']
r = data['r']
rt = data['rt']
Vx = data['Vx']
alpha = data['alpha']
Vf = data['Vf']

U = 2 * np.pi * Rmean * rpm / 60
phi = Vx / U

theta = rt / r
theta = theta[:, :-1, :, :]

na = -1
nb = -2
x1 = x[na, nb, :, :]
t1 = theta[na, nb, :, :]
Vx1 = Vx[na, nb, :, :]
alpha1 = alpha[na]
Vf1 = Vf[nb]

na = 0
# nb = -2
x2 = x[na, nb, :, :]
t2 = theta[na, nb, :, :]
Vx2 = Vx[na, nb, :, :]
alpha2 = alpha[na]
Vf2 = Vf[nb]

fig1, ax1 = plt.subplots(figsize=(10, 5))
fig2, ax2 = plt.subplots(figsize=(10, 5))

levels = np.arange(np.min(Vx2), np.max(Vx1), 5)

cmap = 'seismic'
CS1 = ax1.contourf(t1, x1, Vx1, cmap=cmap, levels=levels)
CS2 = ax2.contourf(t2, x2, Vx2, cmap=cmap, levels=levels)

col = 'r'
ax1.contour(CS1, levels=[0], colors=col)
ax2.contour(CS2, levels=[0], colors=col)

ax1.set_axis_off()
ax2.set_axis_off()

# fig.colorbar(CStop, ax=axtop, location='right')

# fig.subplots_adjust(right=0.8)

# cbar_top = fig.colorbar(CStop)

# ax1.set_title(fr"$ \alpha = {alpha1}^\circ, V = {Vf1} $ m/s")
# ax2.set_title(fr"$ \alpha = {alpha2}^\circ, V = {Vf2} $ m/s")


fig1.set_constrained_layout(True)
fig2.set_constrained_layout(True)

fig1.savefig("Figures/IPM5_CFD_Intake_contours_attached.pdf")
fig2.savefig("Figures/IPM5_CFD_Intake_contours_separated.pdf")
plt.show()