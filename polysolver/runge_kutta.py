"""
Runge-Kutta method for solving
differential equations.


"""

from typing import Callable

def order1(
    fun:Callable,
    x:float,
    y:float,
    z:float,
    h:float
):
    """
    First order Runge-Kutta coefficients.
    
    Use this function to get :math:`k_1` or :math:`l_1`, depending on the function provided.

    Parameters
    ----------
    fun : Callable
        The derivative to use in computing this coefficient.
    x : float
        The x value.
    y : float
        The y value.
    z : float
        The z value.
    h : float
        The step size.
    
    Returns
    -------
    float
        The first order Runge-Kutta coefficient.
    
    Examples
    --------
    >>> k1 = order1(yprime,x,y,z,h)
    >>> l1 = order1(zprime,x,y,z,h)
    """
    return h*fun(x,y,z)
def order2(
    fun:Callable,
    x:float,
    y:float,
    z:float,
    h:float,
    k1:float,
    l1:float
):
    """
    Second order Runge-Kutta coefficients.
    
    Use this function to get :math:`k_2` or :math:`l_2`, depending on the function provided.

    Parameters
    ----------
    fun : Callable
        The derivative to use in computing this coefficient.
    x : float
        The x value.
    y : float
        The y value.
    z : float
        The z value.
    h : float
        The step size.
    k1 : float
        The first order Runge-Kutta coefficient.
    l1 : float
        The first order Runge-Kutta coefficient.
    
    Returns
    -------
    float
        The second order Runge-Kutta coefficient.
    
    Examples
    --------
    >>> k2 = order2(yprime,x,y,z,h,k1,l1)
    >>> l2 = order2(zprime,x,y,z,h,k1,l1)
    """

    return h*fun(
        x+0.5*h,
        y+0.5*k1,
        z+0.5*l1
    )
def order3(
    fun:Callable,
    x:float,
    y:float,
    z:float,
    h:float,
    k2:float,
    l2:float,
):
    """
    Third order Runge-Kutta coefficients.
    
    Use this function to get :math:`k_3` or :math:`l_3`, depending on the function provided.
    
    Parameters
    ----------
    fun : Callable
        The derivative to use in computing this coefficient.
    x : float
        The x value.
    y : float
        The y value.
    z : float
        The z value.
    h : float
        The step size.
    k2 : float
        The second order Runge-Kutta coefficient.
    l2 : float
        The second order Runge-Kutta coefficient.
    
    Returns
    -------
    float
        The third order Runge-Kutta coefficient.
    
    Examples
    --------
    >>> k3 = order3(yprime,x,y,z,h,k2,l2)
    >>> l3 = order3(zprime,x,y,z,h,k2,l2)
    """
    return h*fun(
        x+0.5*h,
        y+0.5*k2,
        z+0.5*l2
    )
def order4(
    fun:Callable,
    x:float,
    y:float,
    z:float,
    h:float,
    k3:float,
    l3:float,
):
    """
    Fourth order Runge-Kutta coefficients.
    
    Use this function to get :math:`k_4` or :math:`l_4`, depending on the function provided.
    
    Parameters
    ----------
    fun : Callable
        The derivative to use in computing this coefficient.
    x : float
        The x value.
    y : float
        The y value.
    z : float
        The z value.
    h : float
        The step size.
    k3 : float
        The third order Runge-Kutta coefficient.
    l3 : float
        The third order Runge-Kutta coefficient.
    
    Returns
    -------
    float
        The fourth order Runge-Kutta coefficient.
    
    Examples
    --------
    >>> k4 = order4(yprime,x,y,z,h,k3,l3)
    >>> l4 = order4(zprime,x,y,z,h,k3,l3)
    """
    return h*fun(
        x+h,
        y+k3,
        z+l3
    )

def dy_and_dz(
    yprime:Callable,
    zprime:Callable,
    x:float,
    y:float,
    z:float,
    h:float
):
    """
    Use the fourth order Runge-Kutta method to compute :math:`\\frac{dy}{dx}` and :math:`\\frac{dz}{dx}`.
    
    Parameters
    ----------
    yprime : Callable
        The :math:`\\frac{dy}{dx}` function.
    zprime : Callable
        The :math:`\\frac{dz}{dx}` function.
    x : float
        The x value.
    y : float
        The y value.
    z : float
        The z value.
    h : float
        The step size.
    
    Returns
    -------
    float, float
        The :math:`\\frac{dy}{dx}` and :math:`\\frac{dz}{dx}` values.
    
    Examples
    --------
    >>> dy,dz = dy_and_dz(yprime,zprime,x,y,z,h)
    """
    k1 = order1(yprime,x,y,z,h)
    l1 = order1(zprime,x,y,z,h)
    k2 = order2(yprime,x,y,z,h,k1,l1)
    l2 = order2(zprime,x,y,z,h,k1,l1)
    k3 = order3(yprime,x,y,z,h,k2,l2)
    l3 = order3(zprime,x,y,z,h,k2,l2)
    k4 = order4(yprime,x,y,z,h,k3,l3)
    l4 = order4(zprime,x,y,z,h,k3,l3)
    return (
        k1/6+k2/3+k3/3+k4/6,
        l1/6+l2/3+l3/3+l4/6
    )

def get_next_xyz(
    yprime:Callable,
    zprime:Callable,
    x:float,
    y:float,
    z:float,
    h:float
):
    """
    Use the fourth order Runge-Kutta method to compute the next x, y, and z values.
    
    Parameters
    ----------
    yprime : Callable
        The :math:`\\frac{dy}{dx}` function.
    zprime : Callable
        The :math:`\\frac{dz}{dx}` function.
    x : float
        The x value.
    y : float
        The y value.
    z : float
        The z value.
    h : float
        The step size.
    
    Returns
    -------
    float, float, float
        The next x, y, and z values.
    
    Examples
    --------
    >>> x,y,z = get_next_xyz(yprime,zprime,x,y,z,h)
    """
    dy,dz = dy_and_dz(yprime,zprime,x,y,z,h)
    dx = h
    return x+dx,y+dy,z+dz