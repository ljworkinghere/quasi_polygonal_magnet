# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 10:17:46 2025

@author: Jie Li
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve
plt.style.use('fivethirtyeight')
mpl.rcParams['text.usetex'] = True

plt.figure(figsize=(10, 4.5))
# range of rho in z-plane
rho_range = np.linspace(0.40, 1, 12)
theta = np.linspace(0, np.pi*2, 1000)

for rho in rho_range:
    x = []
    y = []

    for theta_temp in theta:
        # coordinate in z-plane
        z = rho*np.exp(1j*theta_temp)
        # corresponding coordinate in zeta-plane
        zeta = -z**3 + z + 2/z # the mapping we use
        x.append(zeta.real)
        y.append(zeta.imag)
    plt.subplot(121)
    # the circle in z-plane
    plt.plot(x, y, label=rho, linewidth=1,\
            color = "black")
    plt.subplot(122)
    # the quasi-rectangle in zeta-plane
    plt.plot(np.real(rho*np.exp(1j*theta)),\
              np.imag(rho*np.exp(1j*theta)),\
                linewidth=1, color = "black")

rho  = 0.64 # choose a former we want in z-plane
x = []
y = []
for theta_temp in theta:
    z = rho*np.exp(1j*theta_temp)
    zeta = -z**3 + z + 2/z
    x.append(zeta.real)
    y.append(zeta.imag)
plt.subplot(121)
plt.plot(x, y, label=rho, linewidth=3.0,\
          color = "red")
plt.subplot(122)
plt.plot(np.real(rho*np.exp(1j*theta)),\
          np.imag(rho*np.exp(1j*theta)),\
          linewidth=3.0,\
          color = "red")


plt.subplot(121)
plt.axis("equal")
plt.title(r"$\zeta \left( z \right)$")
# plt.legend()
plt.subplot(122)
mpl.rcParams['text.usetex'] = True
plt.title(r"$z$")
plt.axis("equal")
plt.show()