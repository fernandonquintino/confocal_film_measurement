import numpy as np

from .alpha_correlations import alpha_homo
from .combined_properties import Gw, Jv, rho_alpha, rho_x
from .dimensionless_number import (
    Bond,
    Convection,
    Froude,
    Reynolds,
    Weber,
    martinelli,
)


def friction_f(Re, e=0, d=1.0):
    if Re < 2300:
        f = 64/Re
    elif e==0:
        f = 0.25*(np.log10(150.39/(Re**0.98865) - 152.66/Re))**-2
    else:
        Rr = e / d
        f = 1.613*(np.log(0.234*Rr**1.1007 - 60.525/(Re**1.1105) + 56.291/(Re**1.0712)))**-2
    return f

# Homogeneous based

def dpdz_h(g, d, rho, f):
    dpdz = g**2/(2*d*rho)*f
    return dpdz

def dpdz_mcaddams_1942(g, d, x, rho_l, rho_g, mu_l, mu_g, e=0):
    mu = ((1-x)/mu_l + x/mu_g)**-1
    Re = Reynolds(g, d, mu)
    rho = rho_x(x, rho_l, rho_g)
    f = friction_f(Re, e)
    dpdz = dpdz_h(g, d, rho, f)
    return dpdz

def dpdz_cicchitti_1960(g, d, x, rho_l, rho_g, mu_l, mu_g, e=0):
    mu = (1-x)*mu_l + x*mu_g
    Re = Reynolds(g, d, mu)
    rho = rho_x(x, rho_l, rho_g)
    f = friction_f(Re, e)
    dpdz = dpdz_h(g, d, rho, f)
    return dpdz

def dpdz_dukler_1964(g, d, x, rho_l, rho_g, mu_l, mu_g, e=0):
    alpha_h = alpha_homo(x, rho_l, rho_g)
    mu = (1 - alpha_h)*mu_l + alpha_h*mu_g
    rho = rho_alpha(x, rho_l, rho_g)
    beta = rho_l/rho*(1 - alpha_h) + rho_g/rho*alpha_h
    Re = g*d/mu*beta
    f = friction_f(Re, e)
    dpdz = dpdz_h(g, d, rho, f)
    return dpdz

def dpdz_beatti_walley_1982(g, d, x, rho_l, rho_g, mu_l, mu_g, e=0):
    alpha_h = alpha_homo(x, rho_l, rho_g)
    mu = mu_l*(1-alpha_h)*(1 + 2.5*alpha_h) + mu_g*alpha_h
    Re = Reynolds(g, d, mu)
    rho = rho_x(x, rho_l, rho_g)
    f = friction_f(Re, e)
    dpdz = dpdz_h(g, d, rho, f)
    return dpdz

def dpdz_lin_1991(g, d, x, rho_l, rho_g, mu_l, mu_g, e=0):
    mu = mu_l*mu_g/(mu_g + x**1.4*(mu_l - mu_g))
    Re = Reynolds(g, d, mu)
    rho = rho_x(x, rho_l, rho_g)
    f = friction_f(Re, e)
    dpdz = dpdz_h(g, d, rho, f)
    return dpdz

def dpdz_chen_2001(g, d, x, rho_l, rho_g, mu_l, mu_g, sigma, e=0, ggrav=9.81):
    rho = rho_x(x, rho_l, rho_g)
    Bd = Bond(d, rho_l, rho_g, sigma)
    We = Weber(g, d, rho, sigma)
    if Bd < 10:
        omega = 1 + (0.2 - 0.9*np.exp(-Bd/4))
    else:
        omega = 1 + (We**0.2/(np.exp(Bd/4)**0.3))
    mu = mu_l*mu_g/(mu_g + x**1.4*(mu_l - mu_g))
    Re = Reynolds(g, d, mu)
    f = friction_f(Re, e)
    dpdz = dpdz_h(g, d, rho, f)*omega
    return dpdz

def dpdz_tibirica_2017(g, d, x, rho_l, rho_g, mu_l, mu_g, e=0):
    mu = (1-x)*mu_l + x*mu_g
    alpha_h = alpha_homo(x, rho_l, rho_g)
    rho = rho_l*(1 - alpha_h) + rho_g*alpha_h
    Re = Reynolds(g, d, mu)
    f = 1.415*(rho_l/rho_g)**-0.3263*Re**-0.2342*((rho_l - rho_g)/rho_l)**6.0858
    dpdz = dpdz_h(g, d, rho, f)
    return dpdz

# $\phi$ based

def dpdz_chisholm_1967(g, d, x, rho_l, rho_g, mu_l, mu_g):
    Re_l = Reynolds(g*(1-x), d, mu_l)
    Re_g = Reynolds(g*x, d, mu_g)
    if (Re_l < 2300) & (Re_g < 2300):
        C = 5
    elif (Re_l < 2300) & (Re_g >= 2300):
        C = 12
    elif (Re_l <= 2300) & (Re_g < 2300):
        C = 10
    else:
        C = 20
    X = martinelli(g, d, x, rho_l, rho_g, mu_l, mu_g)
    phi2 = 1 + C/X + 1/(X**2)

    f = friction_f(Re_l)
    dpdz_l = dpdz_h(g*(1-x), d, rho_l, f)
    dpdz = phi2*dpdz_l
    return dpdz

def dpdz_michima_hibiki_1996(g, d, x, rho_l, rho_g, mu_l, mu_g):
    Re_l = Reynolds(g*(1-x), d, mu_l)
    C = 21*(1-np.exp(-0.319*d))
    X = martinelli(g, d, x, rho_l, rho_g, mu_l, mu_g)
    phi2 = 1 + C/X + 1/(X**2)
    f = friction_f(Re_l)
    dpdz_l = dpdz_h(g*(1-x), d, rho_l, f)
    dpdz = phi2*dpdz_l
    return dpdz

def dpdz_sun_mishima_2009(g, d, x, rho_l, rho_g, mu_l, mu_g):
    Re_l = Reynolds(g*(1-x), d, mu_l)
    Re_g = Reynolds(g*x, d, mu_g)
    C = 1.79*(Re_g/Re_l)**0.4*((1-x)/x)**0.5
    X = martinelli(g, d, x, rho_l, rho_g, mu_l, mu_g)
    phi2 = 1 + C/(X**1.9) + 1/(X**2)
    f = friction_f(Re_l)
    dpdz_l = dpdz_h(g*(1-x), d, rho_l, f)
    dpdz = phi2*dpdz_l
    return dpdz

def dpdz_friedel_1979(g, d, x, rho_l, rho_g, mu_l, mu_g, sigma, ggrav=9.78):
    rho = rho_x(x, rho_l, rho_g)
    Fr = Froude(g, d, rho)
    We = Weber(g, d, rho, sigma)
    Re_lo = Reynolds(g, d, mu_l)
    Re_go = Reynolds(g, d, mu_g)
    f_lo = friction_f(Re_lo)
    f_go = friction_f(Re_go)
    E = (1-x)**2 + x**2*(rho_l/rho_g)*(f_go/f_lo)
    F = x**0.78*(1-x)**0.24
    H = (rho_l/rho_g)**0.91*(mu_g/mu_l)**0.19*(1-mu_g/mu_l)**0.70
    phi2 = E + 3.24*F*H/(Fr**0.045*We**0.035)
    f = friction_f(Re_lo)
    dpdz_lo = dpdz_h(g, d, rho_l, f)
    dpdz = phi2*dpdz_lo
    return dpdz

def dpdz_muller_heck_1986(g, d, x, rho_l, rho_g, mu_l, mu_g):
    Re_lo = Reynolds(g, d, mu_l)
    Re_go = Reynolds(g, d, mu_g)
    f_lo = friction_f(Re_lo)
    f_go = friction_f(Re_go)
    dpdz_lo = dpdz_h(g, d, rho_l, f_lo)
    dpdz_go = dpdz_h(g, d, rho_g, f_go)
    Y2 = dpdz_go/dpdz_lo
    phi2 = (1+2*(Y2-1)*x)*(1-x)**(1/3)+Y2*x**3
    dpdz = phi2*dpdz_lo
    return dpdz

def dpdz_nie_2023(g, d, x, rho_l, rho_g, mu_l, mu_g, sigma):
    Bd = Bond(d, rho_l, rho_g, sigma)
    gw = Gw(rho_l, rho_g, d, sigma)
    J = Jv(g, x, d, rho_l, rho_g)
    rho_tp = rho_x(x, rho_l, rho_g)
    Fr = Froude(g, d, rho_tp)
    Re_l = Reynolds(g*(1-x), d, mu_l)
    Re_g = Reynolds(g*x, d, mu_g)
    fl = friction_f(Re_l)
    fg = friction_f(Re_g)
    dpdz_l = dpdz_h(g*(1-x), d, rho_l, fl)
    dpdz_g = dpdz_h(g*x, d, rho_g, fg)
    Co = Convection(x, rho_l, rho_g)
    if (J>2.5) & (g>gw):
        C = 0.94*Fr**0.26*x**-0.4*Bd**0.05*J**-0.35
        n = 0.55
    else:
        C = Co**0.47 + 0.291*Fr**0.55*Bd**0.13*J**-0.74
        n = 0.73
    dpdz = dpdz_g + C*(-dpdz_g**(1-n/2))*(-dpdz_l**(n/2))-dpdz_l
    return dpdz