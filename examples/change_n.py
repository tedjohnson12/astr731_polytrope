from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from polysolver import solve_rust as solve
from polysolver.analysis import central_over_mean_density, xi_1, theta_prime_xi1


if __name__ in '__main__':
    x_init = 0.001
    step_size = 0.001
    
    fig,ax = plt.subplots(1,1,figsize=(10,6))
    for n in np.linspace(0, 4, 5):
        x, y, z = solve(x_init=x_init, n=n, h=step_size, max_iter=100000)
        rho_c_over_rho = central_over_mean_density(np.array(x),np.array(y),n)
        xi1 = xi_1(x)
        theta_prime = theta_prime_xi1(z)

        ax.plot(x, y, label=f'n={n:.2f}\n $\\rho_c / <\\rho> = $ {rho_c_over_rho:.1f}\n $-\\frac{{d\\theta_n}}{{d\\xi}}(\\xi_1) = $ {-theta_prime:.5f} \n $\\xi_1 = ${xi1:.2f}')

    outfile = Path(__file__).parent / 'test_n.png'
    ax.legend()
    ax.set_xlabel('$\\xi$')
    ax.set_ylabel('$\\theta_n$')

    plt.savefig(outfile)
