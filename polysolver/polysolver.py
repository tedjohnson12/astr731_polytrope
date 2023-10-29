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
from typing import Tuple, List
import numpy as np

from polysolver import runge_kutta
from polysolver import derivatives


def solve_python(x_init,n,h,max_iter=1000)->Tuple[List,List]:
    """
    Solve the Lane-Emden equation using a fourth-order Runge-Kutta method.
    
    Parameters
    ----------
    x_init : float
        The initial x value.
    n : int
        The index of the polytrope.
    h : float
        The step size.
    max_iter : int, optional
        The maximum number of iterations. The default is 1000.
    """
    x_prev = x_init
    y_prev = 1
    z_prev = 0
    yprime = derivatives.get_yprime()
    zprime = derivatives.get_zprime(n)
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
    return xs, ys


def solve_rust(
    x_init:float,
    n:int,
    h:float,
    max_iter:int=1000
):
    import polysolver_rust
    return polysolver_rust.solve(x_init,n,h,max_iter)