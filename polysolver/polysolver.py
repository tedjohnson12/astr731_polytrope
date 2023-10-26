"""
Main module for polysolver

The goal of this module is to solve the
Lane-Emden equation:

.. math::
    \\frac{1}{\\xi^2} \\frac{d}{d\\xi}
    \\left( \\xi^2 \\frac{d\\theta_n}{d\\xi} \\right)
    = -\\theta_n^n

Now let's recast this with:

.. math::
    x = \\xi \\\\
    y = \\theta_n \\\\
    z = \\frac{d\\theta_n}{d\\xi} = \\frac{dy}{dx}

We are left with:

.. math::
    y' = \\frac{dy}{dx} = z \\\\
    z' = \\frac{dz}{dx} = -y^n - \\frac{2}{x} z

We will use a fourth-order Runge-Kutta method.

"""

import numpy as np
from typing import Callable

from polysolver import runge_kutta


def yprime(x,y,z):
    return z
def get_zprime(n)->Callable:
    def zprime(x,y,z):
        return -y**n - 2/x*z
    return zprime

def solve(x_init,n,h,max_iter=1000):
    x_prev = x_init
    y_prev = 1
    z_prev = 0
    zprime = get_zprime(n)
    n_iter = 0
    xs = []
    ys = []
    while y_prev > 0 and n_iter < max_iter:
        n_iter += 1
        xs.append(x_prev)
        ys.append(y_prev)
        x_next, y_next, z_next = runge_kutta.get_next_xyz(
            yprime,
            zprime,
            x_prev,
            y_prev,
            z_prev,
            h
        )
        x_prev, y_prev, z_prev = x_next, y_next, z_next
    return np.array(xs), np.array(ys)
    
    