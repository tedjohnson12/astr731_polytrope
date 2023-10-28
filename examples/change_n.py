from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from polysolver import solve


if __name__ in '__main__':
    x_init = 0.01
    step_size = 0.01
    for n in np.linspace(0, 4, 11):
        x, y = solve(x_init=x_init, n=n, h=step_size, max_iter=1000)

        plt.plot(x, y, label=f'n={n:.2f}')

    outfile = Path(__file__).parent / 'test_n.png'

    plt.savefig(outfile)
