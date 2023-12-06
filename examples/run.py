from pathlib import Path
import matplotlib.pyplot as plt
import sys

print(sys.path)
import polysolver
# from polysolver import polysolver_python


if __name__ in '__main__':
    n = 4
    x_init = 0.01
    step_size = 1e-4
    x, y, z = polysolver.solve(x_init=x_init, n=n, h=step_size, max_iter=10000)

    plt.plot(x, y)

    outfile = Path(__file__).parent / 'test_run.png'

    plt.savefig(outfile)
