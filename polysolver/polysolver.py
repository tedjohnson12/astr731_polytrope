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
    zs = []
    while y_prev > 0 and n_iter < max_iter:
        n_iter += 1
        xs.append(x_prev)
        ys.append(y_prev)
        zs.append(z_prev)
        x_next, y_next, z_next = runge_kutta.get_next_xyz(
            yprime,
            zprime,
            x_prev,
            y_prev,
            z_prev,
            h
        )
        x_prev, y_prev, z_prev = x_next, y_next, z_next
    return np.array(xs), np.array(ys), np.array(zs)


def solve_rust(
    x_init:float,
    n:float,
    h:float,
    max_iter:int=1000
):
    """
    Solve the Lane-Emden equation using a fourth-order Runge-Kutta method.
    Implemented in rust.
    
    Parameters
    ----------
    x_init : float
        The initial x value.
    n : float
        The index of the polytrope.
    h : float
        The step size.
    max_iter : int, optional
        The maximum number of iterations. The default is 1000.
    """
    # pylint: disable-next=no-name-in-module
    from polysolver import polysolver_rust
    x,y,z = polysolver_rust.solve(x_init,n,h,max_iter)
    return np.array(x), np.array(y), np.array(z)

def solve(
    x_init:float,
    n:float,
    h:float,
    max_iter:int=1000,
    impl:str='rust'
):
    """
    Solve the Lane-Emden equation using a fourth-order Runge-Kutta method.
    
    Parameters
    ----------
    x_init : float
        The initial x value. Choose something small.
    n : float
        The index of the polytrope.
    h : float
        The step size. This should be less than the pressure scale height.
    max_iter : int, optional
        The maximum number of iterations. The default is 1000.
    impl : str, optional
        The implementation to use. The default is 'rust'.
    
    Returns
    -------
    x : np.ndarray
        The x values. Recall that :math:`x=\\xi`.
    y : np.ndarray
        The y values. Recall that :math:`y=\\theta_n`.
    z : np.ndarray
        The z values. Recall that :math:`z=\\frac{d\\theta_n}{d\\xi}`.
    """
    if impl == 'rust':
        return solve_rust(x_init,n,h,max_iter)
    if impl == 'python':
        return solve_python(x_init,n,h,max_iter)
    else:
        raise NotImplementedError('impl must be "rust" or "python"')