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
X_INIT = 0.1*HLO
N = 1
MAX_ITER = 10000
IMPL = 'rust'

FILENAME = 'res_n1.pdf'
PATH = paths.figures / FILENAME


fig, ax = plt.subplots(2,1,sharex=True,height_ratios=[5,2])
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
    true = Star.analytic(
        x=star.x,
        n=N
    )
    res = (star.y - true.y)/true.y*100
    
    ax[0].plot(star.x,star.y,label=f'{h:.2f}',c=c)
    ax[1].plot(star.x,res,c=c)

ax[1].set_xlabel('$\\xi$')
ax[0].set_ylabel('$\\theta_1$')
ax[1].set_ylabel('Residual (%)')
ax[0].legend()

fig.savefig(PATH)
    
    
    
    