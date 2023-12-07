"""
Resolution Study to find
deviations of :math:`\\rho_c/<\\rho>`
"""

import matplotlib.pyplot as plt
import numpy as np
import paths

from polysolver import Star

plt.style.use('seaborn-v0_8')

HLO = 0.0001
HHI = 0.1
NMODELS = 200
STEPS = np.logspace(
    np.log10(HLO),
    np.log10(HHI),
    NMODELS
)
X_INIT = 1e-20
N = 0
MAX_ITER = 1000000
IMPL = 'rust'

FILENAME = 'res_rho.pdf'
PATH = paths.figures / FILENAME

data = [
    (0., 1.0),
    (1.0, np.pi**2/3),
    (1.5, 5.9907),
    (2.0, 11.402),
    (3.0, 54.183),
    (4.0, 622.41)
]


fig, ax = plt.subplots(1,1,figsize=(5.5,4))
fig.subplots_adjust(right=0.95,top=0.95)


colors = plt.cm.viridis(np.linspace(0,1,len(data)))

for (n, rho), c in zip(data, colors):
    dat = []
    for h in STEPS:
        star = Star.from_soln(
            x_init=X_INIT,
            n=n,
            h=h,
            max_iter=MAX_ITER,
            impl=IMPL
        )
        rho_measured = star.rho_c_over_rho
        dat.append(rho_measured)
    dat = np.array(dat)
    res = np.abs(dat-rho)/rho
    ax.plot(STEPS, res, c=c,label=f'n={n:.1f}')
        

ax.set_xlabel('$h$')
ax.set_xscale('log')
ax.set_ylabel('$\\Delta \\frac{\\rho_c}{<\\rho>}$')
ax.set_yscale('log')
ax.legend()

fig.savefig(PATH)
    
    
    
    