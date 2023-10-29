from pathlib import Path
import matplotlib.pyplot as plt
import sys
from datetime import datetime

print(sys.path)
import polysolver

if __name__ in '__main__':
    fig,ax = plt.subplots(1,2,figsize=(10,5))
    langs = ['Python','Rust']
    for i, solve in enumerate([polysolver.solve_python, polysolver.solve_rust]):
        n = 4
        x_init = 1e-6
        step_size = 1e-6
        start_time = datetime.now()
        x, y = solve(x_init=x_init, n=n, h=step_size, max_iter=10000000)
        dtime = datetime.now() - start_time
        ax[i].plot(x, y)
        ax[i].set_title(f'{langs[i]} $\\Delta t = ${dtime.total_seconds():.3f} s')
    
    outfile = Path(__file__).parent / 'test_langs.png'

    plt.savefig(outfile)