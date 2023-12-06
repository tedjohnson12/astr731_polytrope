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

def get_rho_norm(y:np.ndarray,n:float):
    rho_over_rho_c = y**n
    return rho_over_rho_c

def norm_mass(x:np.ndarray,y:np.ndarray,n:float):
    r = x
    dv = 4 * np.pi * r**2
    rho_over_rho_c = get_rho_norm(y,n)
    dm_over_rho_c = dv*rho_over_rho_c
    mass_over_rho_c = np.trapz(dm_over_rho_c,r)
    return mass_over_rho_c
def volume(x:np.ndarray):
    return 4/3 *np.pi * xi_1(x)**3
    
def xi_1(x):
    return x[-1]

def theta_prime_xi1(z):
    return z[-1]

def central_over_mean_density(x,y,n):
    """
    .. math::
        r = r_n x \\\\
        \\rho = \\rho_c y^n
        
    """
    mass_over_rhoc = norm_mass(x,y,n)
    vol = volume(x)
    return vol/mass_over_rhoc