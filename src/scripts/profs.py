"""
Make a plot of profiles for different values of n.
"""

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
import paths

from polysolver import Star

plt.style.use('seaborn-v0_8')

H = 1e-3
NMODELS = 1000
NS = np.linspace(0,4,NMODELS)
X_INIT = 1e-20
MAX_ITER = 100000
IMPL = 'rust'
N_RESAMPLE = 100
ALPHA = 0.1
FILENAME = 'profs.pdf'
PATH = paths.figures / FILENAME

fig, ax = plt.subplots(1,1,figsize=(5.5,4))
fig.subplots_adjust(right=0.95,top=0.95)

colors = plt.cm.viridis(np.linspace(0,1,NMODELS))
mappable = ScalarMappable(
    norm=Normalize(0,4),
    cmap=plt.cm.viridis
)


for n,c in zip(NS,colors):
    star = Star.from_soln(
        x_init=X_INIT,
        n=n,
        h=H,
        max_iter=MAX_ITER,
        impl=IMPL
    )
    xi1 = star.xi1
    xnew = np.linspace(X_INIT, xi1, N_RESAMPLE)
    ynew = star.resample_y(xnew)
    pcen = star.central_pressure(mass=1, radius=1)
    pressure = ynew**(1+n)*pcen
    ax.plot(xnew/xi1, pressure, label=f'n={n:.2f}',c=c,alpha=ALPHA)

ax.set_xlabel('$r/R$')
ax.set_ylabel('Pressure (dyne cm$^{-2}$)')
ax.set_yscale('log')

fig.colorbar(mappable, ax=ax, label='n')

fig.savefig(PATH)