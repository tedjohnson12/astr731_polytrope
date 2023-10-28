"""

Derivates to use in the Runge-Kutta method.

"""
from typing import Callable


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
    def zprime(x, y, z):
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
        return -y**n - 2/x*z
    return zprime
