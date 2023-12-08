"""
Recreate table 7.1

"""
import numpy as np

from polysolver import Star
import paths



data = {
    0.0:
        {
            'xi1': np.sqrt(6),
            'dtheta': np.sqrt(6)/3,
            'rho': 1
        },
    1.0:
        {
            'xi1': np.pi,
            'dtheta': 1/np.pi,
            'rho': np.pi**2/3
        },
    1.5:
        {
            'xi1': 3.6538,
            'dtheta': 0.20330,
            'rho': 5.9907
        },
    2.0:
        {
            'xi1': 4.3529,
            'dtheta': 0.12725,
            'rho': 11.402
        },
    3.0:
        {
            'xi1':6.8969,
            'dtheta': 0.04243,
            'rho': 54.183
        },
    4.0:
        {
            'xi1': 14.972,
            'dtheta': 0.00802,
            'rho': 662.41
        }
}

NMODELS = 4*4 + 1
NS = np.linspace(0,4,NMODELS)
H = 1e-3
X_INIT = 1e-20

NAMES = [
    '$n$',
    '$\\xi_1$',
    'exp.',
    '$-\\theta_n\\\'(\\xi_1)$',
    'exp.',
    '$\\rho_c / <\\rho>$',
    'exp.'
]
TITLE = ' & '.join(NAMES)

filename = 'tab7.txt'
outfile = paths.output / filename

CAPTION = """\
Relevent quantities for various polytropic
solutions to the Lane-Emden equation. The 
column exp. lists the value for the previous column
given by \\citet{textbook}.
"""

def format(x:float,n:int):
    return f'{x:.{n}g}'


def line(
    n:float,
    xi1:float,
    dtheta:float,
    rho:float
):
    entries = [
        format(n,3),
        format(xi1,5),
        format(data[n]['xi1'],5) if n in data else '-',
        format(dtheta,5),
        format(data[n]['dtheta'],5) if n in data else '-',
        format(rho,5),
        format(data[n]['rho'],5) if n in data else '-',
    ]
    return ' & '.join(entries)

with open(outfile,'w',encoding='utf-8') as f:
    f.write(r'\begin{table}[ht]'+'\n')
    f.write(r'\centering'+'\n')
    f.write(r'\begin{tabular}{ccccccc}'+'\n')
    f.write(r'\hline'+'\n')
    f.write(f'{TITLE} \\\\'+'\n')
    f.write(r'\hline' + '\n')
    for n in NS:
        star = Star.from_soln(
            x_init=X_INIT,
            n=n,
            h=H,
            max_iter=100000,
            impl='rust'
        )
        f.write(
            line(
                n,
                star.xi1,
                star.theta_prime,
                star.rho_c_over_rho
            ) + ' \\\\' + '\n'
        )
    f.write(r'\hline' +'\n')
    f.write(r'\end{tabular}' + '\n')
    f.write(r'\caption{' + CAPTION + '}' + '\n')
    f.write(r'\label{tab:7}' + '\n')
    f.write(r'\end{table}' + '\n')
    