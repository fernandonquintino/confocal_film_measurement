def Reynolds(g, d, mu):
    Re = g*d/mu
    return Re

def Convection(x, rho_l, rho_g):
    Co = (1/x-1)**0.8*(rho_g/rho_l)**0.5
    return Co

def Froude(g, d, rho, ggrav=9.78):
    Fr = g**2/(rho**2*ggrav*d)
    return Fr

def Bond(d, rho_l, rho_g, sigma, ggrav=9.78):
    Bd = ggrav*(rho_l - rho_g)*d**2/sigma
    return Bd

def Weber(g: float, d: float, rho: float, sigma: float):
    We = g**2*d/(rho*sigma)
    return We

def martinelli(g, d, x, rho_l, rho_g, mu_l, mu_g):
    Re_l = Reynolds(g*(1-x), d, mu_l)
    Re_g = Reynolds(g*x, d, mu_g)
    if (Re_l < 2300) & (Re_g < 2300):
        X = ((1-x)/x)**0.5*(rho_g/rho_l)**0.5*(mu_l/mu_g)**0.5
    elif (Re_l < 2300) & (Re_g >= 2300):
        X = 18.65*Re_g**-0.4*((1-x)/x)**0.5*(rho_g/rho_l)**0.5*(mu_l/mu_g)**0.5
    elif (Re_l <= 2300) & (Re_g < 2300):
        X = 0.0536*Re_l**0.4*((1-x)/x)**0.5*(rho_g/rho_l)**0.5*(mu_l/mu_g)**0.5
    else:
        X = ((1-x)/x)**0.9*(rho_g/rho_l)**0.5*(mu_l/mu_g)**0.1
    return X