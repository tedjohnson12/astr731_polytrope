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
    dy,dz = dy_and_dz(yprime,zprime,x,y,z,h)
    dx = h
    return x+dx,y+dy,z+dz