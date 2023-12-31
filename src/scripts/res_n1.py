"""
Resolution Study for n=0
"""

import matplotlib.pyplot as plt
import numpy as np
import paths

from polysolver import Star

plt.style.use('seaborn-v0_8')

HLO = 0.01
HHI = 0.1
NMODELS = 5
STEPS = np.logspace(
    np.log10(HLO),
    np.log10(HHI),
    NMODELS
)
X_INIT = 1e-20
N = 1
MAX_ITER = 10000
IMPL = 'rust'
NPOINTS_RESAMPLE = 100

FILENAME = 'res_n1.pdf'
PATH = paths.figures / FILENAME


fig, ax = plt.subplots(2,1,sharex=True,height_ratios=[5,3],figsize=(5.5,4))
fig.subplots_adjust(hspace=0.05,right=0.95,top=0.95)

colors = plt.cm.viridis(np.linspace(0,1,NMODELS))

for h, c in zip(STEPS, colors):
    star = Star.from_soln(
        x_init=X_INIT,
        n=N,
        h=h,
        max_iter=MAX_ITER,
        impl=IMPL
    )
    xi1 = star.xi1
    xnew = np.linspace(X_INIT,xi1,NPOINTS_RESAMPLE)
    ynew = star.resample_y(xnew)
    true = Star.analytic(
        x=xnew,
        n=N
    )
    res = np.abs(ynew - true.y)
    
    ax[0].plot(xnew,ynew,label=f'{h:.2f}',c=c)
    ax[1].plot(xnew,res,c=c)

ax[1].set_xlabel('$\\xi$')
ax[0].set_ylabel('$\\theta_1$')
ax[1].set_ylabel('$|\\theta_1 - \\theta_{1,\\mathrm{true}}|$')
ax[1].set_yscale('log')
ax[0].legend()

fig.savefig(PATH)
    
    
    
    