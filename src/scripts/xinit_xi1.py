"""
Resolution Study to find
deviations of \\xi_1
"""

import matplotlib.pyplot as plt
import numpy as np
import paths

from polysolver import Star

plt.style.use('seaborn-v0_8')

XLO = 1e-10
XHI = 0.1
NMODELS = 200
X_INITS = np.logspace(
    np.log10(XLO),
    np.log10(XHI),
    NMODELS
)
H = 1e-3
N = 0
MAX_ITER = 1000000
IMPL = 'rust'

FILENAME = 'xinit_xi1.pdf'
PATH = paths.figures / FILENAME

data = [
    (0., np.sqrt(6)),
    (1.0, np.pi),
    (1.5, 3.6538),
    (2.0, 4.3529),
    (3.0, 6.8969),
    (4.0, 14.972)
]


fig, ax = plt.subplots(1,1,figsize=(5.5,4))
fig.subplots_adjust(right=0.95,top=0.95)


colors = plt.cm.viridis(np.linspace(0,1,len(data)))

for (n, xi1), c in zip(data, colors):
    dat = []
    for x_init in X_INITS:
        star = Star.from_soln(
            x_init=x_init,
            n=n,
            h=H,
            max_iter=MAX_ITER,
            impl=IMPL
        )
        xi1_measured = star.xi1
        dat.append(xi1_measured)
    dat = np.array(dat)
    res = np.abs(xi1 - dat)/xi1
    ax.plot(X_INITS, res, c=c,label=f'n={n:.1f}')
        

ax.set_xlabel('x_init')
ax.set_xscale('log')
ax.set_ylabel('$\\Delta \\xi_1$')
ax.set_yscale('log')
ax.legend()

fig.savefig(PATH)
    
    
    
    