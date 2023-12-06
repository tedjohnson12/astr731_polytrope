"""
Analysis module for this polytrope solver.

Let's first recall some math.

.. math::
    x = \\xi = \\frac{r}{r_n}

where

.. math::
    r_n^2 = \\frac{(n+1) P_c}{4\\pi G \\rho_c^2}

is the square of the scale length.

We also have

.. math::
    \\rho(r) = \\rho_c \\theta^n(r)

where
.. math::
    \\theta = y


"""
import numpy as np

def get_rho_norm(y:np.ndarray,n:float)->np.ndarray:
    """
    Get the density as a fraction of the central density.
    
    Parameters
    ----------
    y : np.ndarray
        The y values. Recall that :math:`y=\\theta_n`
        and :math:`\\rho(r) = \\rho_c \\theta^n(r)`
    n : float
        The index of the polytrope.
    
    Returns
    -------
    np.ndarray
        The density as a fraction of the central density.
    """
    rho_over_rho_c = y**n
    return rho_over_rho_c

def norm_mass(x:np.ndarray,y:np.ndarray,n:float)->float:
    """
    Get the mass divided by the central density.
    
    Parameters
    ----------
    x : np.ndarray
        The x values. Recall that :math:`x=\\xi=\\frac{r}{r_n}`
    y : np.ndarray
        The y values. Recall that :math:`y=\\theta_n`
        and :math:`\\rho(r) = \\rho_c \\theta^n(r)`
    n : float
        The index of the polytrope.
    
    Returns
    -------
    float
        The mass divided by the central density.
    """
    r = x
    dv = 4 * np.pi * r**2
    rho_over_rho_c = get_rho_norm(y,n)
    dm_over_rho_c = dv*rho_over_rho_c
    mass_over_rho_c = np.trapz(dm_over_rho_c,r)
    return mass_over_rho_c
def volume(x:np.ndarray)->float:
    """
    Get the volume in units of x.
    
    Parameters
    ----------
    x : np.ndarray
        The x values. Recall that :math:`x=\\xi=\\frac{r}{r_n}`
    
    Returns
    -------
    float
        The volume in units of x.
    """
    return 4/3 *np.pi * xi_1(x)**3
    
def xi_1(x):
    """
    Get the last value :math:`\\xi_1`.
    
    Parameters
    ----------
    x : np.ndarray
        The x values. Recall that :math:`x=\\xi=\\frac{r}{r_n}`
    
    Returns
    -------
    float
        :math:`\\xi_1`.
    """
    return x[-1]

def theta_prime_xi1(z):
    """
    Get :math:`-\\frac{d\\theta_n}{d\\xi}(\\xi_1)`.
    """
    return z[-1]

def central_over_mean_density(x,y,n):
    """
    Get the central density as a fraction of the mean density.
    
    Parameters
    ----------
    x : np.ndarray
        The x values. Recall that :math:`x=\\xi=\\frac{r}{r_n}`
    y : np.ndarray
        The y values. Recall that :math:`y=\\theta_n`
        and :math:`\\rho(r) = \\rho_c \\theta^n(r)`
    n : float
        The index of the polytrope.
    
    Returns
    -------
    float
        The central density as a fraction of the mean density.
    
    Notes
    -----
    .. math::
        r = r_n x \\\\
        \\rho = \\rho_c y^n
        
    """
    mass_over_rhoc = norm_mass(x,y,n)
    vol = volume(x)
    return vol/mass_over_rhoc