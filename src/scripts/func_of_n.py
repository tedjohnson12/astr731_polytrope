"""
How the quantities change as 
a function of n
"""

import numpy as np
import matplotlib.pyplot as plt

from polysolver import Star
import paths


H = 1e-3
X_INIT = 1e-20
NMODELS = 100

NS = np.linspace(0,4,NMODELS)

filename = 'func_of_n.pdf'
path = paths.figures / filename

fig, ax = plt.subplots(1,1,figsize=(5.5,4))
fig.subplots_adjust(right=0.95,top=0.95)

xis = np.zeros_like(NS,dtype=np.float64)
dthetas = np.zeros_like(NS,dtype=np.float64)
rhos = np.zeros_like(NS,dtype=np.float64)

for i, n in enumerate(NS):
    star = Star.from_soln(
        x_init=X_INIT,
        n=n,
        h=H,
        max_iter=100000,
        impl='rust'
    )
    xis[i] = star.xi1
    dthetas[i] = star.theta_prime
    rhos[i] = star.rho_c_over_rho

ax.plot(NS,xis,label=r'$\xi_1$')
ax.plot(NS,dthetas,label='$-\\frac{d\\theta_n}{d\\xi}(\\xi_1)$')
ax.plot(NS,rhos,label=r'$\rho_c / <\rho>$')

ax.set_xlabel('n')
ax.set_yscale('log')
ax.legend()

fig.savefig(path)
