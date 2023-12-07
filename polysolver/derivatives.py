"""

Derivates to use in the Runge-Kutta method.

"""
from typing import Callable
import cmath


def get_yprime() -> Callable:
    """
    Get the yprime function.

    Returns
    -------
    Callable
        The yprime function
    """
    def yprime(x, y, z):
        """
        :math:`\\frac{dy}{dx}`

        Parameters
        ----------
        x : float
            The x value.
        y : float
            The y value.
        z : float
            The z value.

        Returns
        -------
        float
            The yprime value.
        """
        return z
    return yprime


def get_zprime(n) -> Callable:
    """
    Get the zprime function.

    Parameters
    ----------
    n : int
        The index of the polytrope.

    Returns
    -------
    Callable
        The zprime function
    """
    def zprime(
        x:float,
        y:float,
        z:float
    ):
        """
        :math:`\\frac{dz}{dx}`

        Parameters
        ----------
        x : float
            The x value.
        y : float
            The y value.
        z : float
            The z value.

        Returns
        -------
        float
            The zprime value.
        """
        if isinstance(y, complex):
            raise RuntimeError('Complex y passed to zprime')
        a = y**n
        if isinstance(a, complex):
            a = -abs(a)
        return -a - 2/x*z
    return zprime
