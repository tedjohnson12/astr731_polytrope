"""
Resolution Study to find
deviations of \\xi_1
"""

import matplotlib.pyplot as plt
import numpy as np
import paths

from polysolver import Star

plt.style.use('seaborn-v0_8')

HLO = 0.0001
HHI = 0.1
NMODELS = 40
STEPS = np.logspace(
    np.log10(HLO),
    np.log10(HHI),
    NMODELS
)
X_INIT = 1e-8
N = 0
MAX_ITER = 1000000
IMPL = 'rust'

FILENAME = 'res_xi1.pdf'
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
    for h in STEPS:
        star = Star.from_soln(
            x_init=X_INIT,
            n=n,
            h=h,
            max_iter=MAX_ITER,
            impl=IMPL
        )
        xi1_measured = star.xi1
        dat.append(xi1_measured)
    dat = np.array(dat)
    res = (xi1 - dat)/xi1
    ax.plot(STEPS, res, c=c,label=f'n={n:.1f}')
        

ax.set_xlabel('$h$')
ax.set_xscale('log')
ax.set_ylabel('$\\frac{\\xi_1 - \\xi_{1, \\text{exp}}}{\\xi_{1, \\text{exp}}}$')
ax.set_yscale('log')
ax.legend()

fig.savefig(PATH)
    
    
    
    