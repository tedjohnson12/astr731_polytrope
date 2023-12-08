"""
Resolution Study to find
deviations of \\xi_1
"""

import matplotlib.pyplot as plt
import numpy as np
import paths

from polysolver import Star

plt.style.use('seaborn-v0_8')

XLO = 1e-8
XHI = 0.1
NMODELS_X = 40
X_INITS = np.logspace(
    np.log10(XLO),
    np.log10(XHI),
    NMODELS_X
)
HLO = 0.0001
HHI = 0.1
NMODELS_H = 30
STEPS = np.logspace(
    np.log10(HLO),
    np.log10(HHI),
    NMODELS_H
)
N = 0
MAX_ITER = 1000000
IMPL = 'rust'

FILENAME = 'xinit_res_xi1.pdf'
PATH = paths.figures / FILENAME

data = [
    (0., np.sqrt(6)),
    (1.0, np.pi),
    (1.5, 3.6538),
    (2.0, 4.3529),
    (3.0, 6.8969),
    (4.0, 14.972)
]


fig, axes = plt.subplots(2,3,figsize=(8,5))
axes = np.ravel(axes)
fig.subplots_adjust(right=0.95,top=0.95,hspace=0.4,wspace=0.4)


colors = plt.cm.viridis(np.linspace(0,1,len(data)))



for (n, xi1), c, ax in zip(data, colors, axes):
    dat = np.zeros((NMODELS_X, NMODELS_H))
    for i,x_init in enumerate(X_INITS):
        for j,h in enumerate(STEPS):
            star = Star.from_soln(
                x_init=x_init,
                n=n,
                h=h,
                max_iter=MAX_ITER,
                impl=IMPL
            )
            xi1_measured = star.xi1
            dat[i,j] = xi1_measured
    res = np.abs(xi1 - dat)/xi1
    
    ax.set_title(f'$\\Delta \\xi_1$ for n={n:.1f}')
    im = ax.pcolormesh(X_INITS, STEPS, np.log10(res.T), cmap='viridis')
    fig.colorbar(im, ax=ax, pad=0.01)
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.axvline(1.0,ls='--',c='k')
    ax.plot(STEPS,STEPS,c='k',ls='--')

fig.text(0.5, 0.04, 'x_init', ha='center')
fig.text(0.04, 0.5, '$h$', va='center', rotation='vertical')


fig.savefig(PATH)
    
    
    
    