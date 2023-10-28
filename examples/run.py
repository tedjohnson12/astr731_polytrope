from pathlib import Path
import matplotlib.pyplot as plt

from polysolver import solve


if __name__ in '__main__':
    n = 0
    x_init = 0.01
    step_size = 0.01
    x, y = solve(x_init=x_init, n=n, h=step_size, max_iter=1000)

    plt.plot(x, y)

    outfile = Path(__file__).parent / 'test_run.png'

    plt.savefig(outfile)