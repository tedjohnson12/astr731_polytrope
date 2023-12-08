"""
x_init Study for n=0
"""

import matplotlib.pyplot as plt
import numpy as np
import paths

from polysolver import Star

plt.style.use('seaborn-v0_8')

XLO = 1e-10
XHI = 0.1
NMODELS = 10
X_INITS = np.logspace(
    np.log10(XLO),
    np.log10(XHI),
    NMODELS
)
H = 1e-3
N = 0
MAX_ITER = 10000
IMPL = 'rust'
NPOINTS_RESAMPLE = 100

FILENAME = 'xinit_n0.pdf'
PATH = paths.figures / FILENAME


fig, ax = plt.subplots(2,1,sharex=True,height_ratios=[5,2],figsize=(5.5,4))
fig.subplots_adjust(hspace=0.05,right=0.95,top=0.95)


colors = plt.cm.viridis(np.linspace(0,1,NMODELS))

for x_init, c in zip(X_INITS, colors):
    star = Star.from_soln(
        x_init=x_init,
        n=N,
        h=H,
        max_iter=MAX_ITER,
        impl=IMPL
    )
    xi1 = star.xi1
    xnew = np.linspace(x_init,xi1,NPOINTS_RESAMPLE)
    ynew = star.resample_y(xnew)
    true = Star.analytic(
        x=xnew,
        n=N
    )
    res = np.abs(ynew - true.y)
    
    ax[0].plot(xnew,ynew,label=f'{x_init:.2e}',c=c)
    ax[1].plot(xnew,res,c=c)

ax[1].set_xlabel('$\\xi$')
ax[0].set_ylabel('$\\theta_0$')
ax[1].set_ylabel('$|\\theta_0 - \\theta_{0,\\mathrm{true}}|$')
ax[1].set_yscale('log')
ax[0].legend()

fig.savefig(PATH)
    
    
    
    