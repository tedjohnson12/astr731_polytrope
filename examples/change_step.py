from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from polysolver import solve_rust as solve


if __name__ in '__main__':
    n = 4
    x_init = 0.01
    for step_size in np.linspace(0.01, 0.4, 11):
        x, y, _ = solve(x_init=x_init, n=n, h=step_size, max_iter=1000)

        plt.plot(x, y, label=f'log($\\Delta \\xi$)={np.log10(step_size):.2f}')

    outfile = Path(__file__).parent / 'test_step.png'

    plt.savefig(outfile)
