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
import warnings
import numpy as np
from scipy.interpolate import interp1d, CubicSpline

from polysolver import solve

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
    xi1 = xi_1(x,y)
    x = np.append(x[:-1],xi1)
    y = np.append(y[:-1],0)
    r = x
    dv = 4 * np.pi * r**2
    rho_over_rho_c = get_rho_norm(y,n)
    dm_over_rho_c = dv*rho_over_rho_c
    mass_over_rho_c = np.trapz(dm_over_rho_c,r)
    return mass_over_rho_c
def volume(x:np.ndarray,y:np.ndarray)->float:
    """
    Get the volume in units of x.
    
    Parameters
    ----------
    x : np.ndarray
        The x values. Recall that :math:`x=\\xi=\\frac{r}{r_n}`
    y : np.ndarray
        The y values. Recall that :math:`y=\\theta_n`
        and :math:`\\rho(r) = \\rho_c \\theta^n(r)`
    
    Returns
    -------
    float
        The volume in units of x.
    """
    return 4/3 *np.pi * xi_1(x,y)**3
    
def xi_1(x,y):
    """
    Get the last value :math:`\\xi_1`.
    
    Parameters
    ----------
    x : np.ndarray
        The x values. Recall that :math:`x=\\xi=\\frac{r}{r_n}`
    y : np.ndarray
        The y values. Recall that :math:`y=\\theta_n`
        and :math:`\\rho(r) = \\rho_c \\theta^n(r)`
    
    Returns
    -------
    float
        :math:`\\xi_1`.
    
    Notes
    -----
    We assume here that the last y value is less than 0.
    Then we use linear interpolation to find the value
    of x at which y=0.
    """
    interp = CubicSpline(-y[-3:],x[-3:])
    return interp(0)

def theta_prime_xi1(z,y):
    """
    Get :math:`-\\frac{d\\theta_n}{d\\xi}(\\xi_1)`.
    """
    interp = CubicSpline(-y[-3:],z[-3:])
    return -interp(0)

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
    vol = volume(x,y)
    return vol/mass_over_rhoc

class Star:
    """
    A polytropic star
    
    Parameters
    ----------
    x : np.ndarray
        The x values. Recall that :math:`x=\\xi=\\frac{r}{r_n}`
    y : np.ndarray
        The y values. Recall that :math:`y=\\theta_n`
        and :math:`\\rho(r) = \\rho_c \\theta^n(r)`
    z : np.ndarray
        The z values. Recall that :math:`z=\\frac{d\\theta_n}{d\\xi}`
    n : float
        The index of the polytrope.
    """
    def __init__(
        self,
        x:np.ndarray,
        y:np.ndarray,
        z:np.ndarray,
        n:float
    ):
        self.x = x
        self.y = y
        self.z = z
        self.n = n
    @classmethod
    def from_soln(
        cls,
        x_init:float,
        n:float,
        h:float,
        max_iter:int=1000,
        impl:str='rust'
    ):
        """
        Create a star from a solution to the Lane-Emden equation.
        
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
        Star
            The star.
        """
        x,y,z = solve(x_init,n,h,max_iter,impl)
        return cls(x,y,z,n)
    @classmethod
    def _zero(cls,x:np.ndarray):
        """
        Analytic solution to the Lane-Emden equation
        for :math:`n=0`.
        
        Parameters
        ----------
        x : np.ndarray
            The x values. Recall that :math:`x=\\xi=\\frac{r}{r_n}`
        
        Notes
        -----
        This is the constant density model.
        
        .. math::
            y = 1 - \\frac{x^2}{6}
        
        .. math::
            z = \\frac{-x}{3}
        
        .. math::
            \\xi_1 = \\sqrt{6}
        """
        xmin = 0
        xmax = np.sqrt(6)
        if np.any(x<xmin):
            warnings.warn(f'Analytic solution only valid for x>{xmin:.2f}',RuntimeWarning)
        if np.any(x>xmax):
            warnings.warn(f'Analytic solution only valid for x<={xmax:.2f}',RuntimeWarning)
        # x = np.where((x<xmin)|(x>xmax),np.nan,x)
        y = 1 - x**2/6
        z = -x/3
        return cls(x,y,z,0)
    @classmethod
    def _one(cls,x:np.ndarray):
        """
        Analytic solution to the Lane-Emden equation
        for :math:`n=1`.
        
        Parameters
        ----------
        x : np.ndarray
            The x values. Recall that :math:`x=\\xi=\\frac{r}{r_n}`
        
        Notes
        -----
        This is the sinc function.
        Note that the numpy implementation
        of sinc implicitly multiplies by :math:`\\pi`.
        
        .. math::
            y = \\frac{\\sin(x)}{x}
        
        .. math::
            z = \\frac{\\cos(x)-\\sin(x)/x}{x}
        
        except :math:`z=0` for :math:`x=0`.
            
        .. math::
            \\xi_1 = \\pi
        """
        xmin = 0
        xmax = np.pi
        if np.any(x<xmin):
            warnings.warn(f'Analytic solution only valid for x>{xmin:.2f}',RuntimeWarning)
        if np.any(x>xmax):
            warnings.warn(f'Analytic solution only valid for x<={xmax:.2f}',RuntimeWarning)
        # x = np.where((x<xmin)|(x>xmax),np.nan,x)
        y = np.sinc(x/np.pi)
        z = np.where(x==0,0,(np.cos(x)-np.sinc(x/np.pi))/x)
        return cls(x,y,z,1)
    @classmethod
    def analytic(
        cls,
        x:np.ndarray,
        n:int
    ):
        """
        Get the analytic solution to the Lane-Emden equation.
        
        Parameters
        ----------
        x : np.ndarray
            The x values. Recall that :math:`x=\\xi=\\frac{r}{r_n}`
        n : float
            The index of the polytrope.
        
        Returns
        -------
        Star
            The star.
        """
        if n==0:
            return cls._zero(x)
        elif n==1:
            return cls._one(x)
        else:
            raise NotImplementedError(f'There is no analytic solution for n={n:d}')
    @property
    def xi1(self)->float:
        """
        Get the value of :math:`\\xi` at the surface.
        """
        return xi_1(self.x,self.y)
    @property
    def theta_prime(self)->float:
        """
        Get :math:`-\\frac{d\\theta_n}{d\\xi}` at the surface.
        """
        return theta_prime_xi1(self.z,self.y)
    @property
    def rho_c_over_rho(self)->float:
        """
        Get the central density divided by the
        mean density.
        """
        return central_over_mean_density(self.x,self.y,self.n)
    def resample_y(
        self,
        x:np.ndarray
    ):
        """
        Resample the y values.
        """
        interp = CubicSpline(self.x,self.y)
        return interp(x)
    def central_pressure(
        self,
        mass:float,
        radius:float
    ):
        """
        Get the central pressure in
        dyne cm-2
        
        Parameters
        ----------
        mass : float
            The mass of the star in solar masses.
        radius : float
            The radius of the star in solar radii.
        
        Returns
        -------
        float
            The central pressure in dyne cm-2.
        """
        num = 8.952e14
        den = (self.n+1)*self.theta_prime**2
        return num/den * mass**2 * radius**-4