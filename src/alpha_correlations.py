import pandas as pd

from .dimensionless_number import Reynolds, Weber


def alpha_homo(x, rho_l, rho_g):
    alpha = (1+(1-x)/x*(rho_g/rho_l))**-1
    return alpha

def alpha_lockhart_martinelli_1949(x, rho_l, rho_g, mu_l, mu_g):
    Sl = 0.28*((1-x)/x)**-0.361*(rho_g/rho_l)**-0.645*(mu_l/mu_g)**0.071
    alpha = (1+(1-x)/x*(rho_g/rho_l)*Sl)**-1
    return alpha

def alpha_zivi_1964(x, rho_l, rho_g):
    Sl = (rho_g/rho_l)**(-1/3)
    alpha = (1+(1-x)/x*(rho_g/rho_l)*Sl)**-1
    return alpha

def alpha_premoli_1970(x, rho_l, rho_g, mu_l, sigma, g, d):
    We_lo = Weber(g, d, rho_l, sigma)
    Re_lo = Reynolds(g, d, mu_l)
    E1 = 1.578*Re_lo**-0.19*(rho_l/rho_g)**0.22
    E2 = 0.0273*We_lo*Re_lo**-0.51*(rho_l/rho_g)**-0.08
    alpha_h = alpha_homo(x, rho_l, rho_g)
    y = alpha_h/(1-alpha_h)
    Sl = 1+E1*((y/(1+y*E2))-y*E2)**0.5
    alpha = (1+(1-x)/x*(rho_g/rho_l)*Sl)**-1
    return alpha

def alpha_kanizawa_ribastki_2015(x, rho_l, rho_g, mu_l, mu_g, g, d, grav=9.81):
    Fr_dp = g**2/((rho_l-rho_g)**2*grav*d)
    Sl = 1.021*Fr_dp**-0.092*((1-x)/x)**(-1/3)*(rho_g/rho_l)**(-2/3)*(mu_l/mu_g)**-0.368
    alpha = (1+(1-x)/x*(rho_g/rho_l)*Sl)**-1
    return alpha

def alpha_tibirica_2017(x, rho_l, rho_g, g, d, grav=9.81):
    Fr_dp = g**2*((rho_l-rho_g)**2*grav*d)**-1
    Sl = 1.2364*Fr_dp**-0.1082*(rho_g/rho_l)**-0.31*((1-x)/x)**-0.267
    alpha = (1+(1-x)/x*(rho_g/rho_l)*Sl)**-1
    return alpha

def alpha_hughmark_1962(x, rho_l, rho_g, mu_l, mu_g, d, g, grav=9.81, max_iterations=100, tol=1e-5):
    alpha_h = alpha_homo(x, rho_l, rho_g)
    alpha = alpha_h
    for i in range(max_iterations):
        Z = (g*d / (mu_l*(1-alpha) + mu_g*alpha) )**(1/6) * ((g**2*x**2) / (grav*d*rho_g**2*alpha_h**2*(1-alpha_h)**2))**(1/8)
        if Z<10:
            K = -0.164 + 0.310*Z - 0.0353*Z**2 + 0.00137*Z**3
        else:
            K = 0.755 + 0.00358*Z - 0.144*10**-4*Z**2
        alpha_new = K*alpha_h

        if abs(alpha_new-alpha)<tol:
            return alpha_new
        alpha = alpha_new 
    print("Warning: Maximum iterations reached. Result may not have converged.")
    return alpha

def alpha_rouhani_1970(x, rho_l, rho_g, mu_l, mu_g, sigma, d, g, grav=9.81):
    C0 = 1.1
    u_gj = 1.18*(1-x)*(grav*sigma*(rho_l-rho_g)/(rho_l**2))**0.25
    alpha = (x/rho_g)/(C0*((1-x)/rho_l+x/rho_g)+(u_gj/g))
    return alpha

def alpha_steiner_1993(x, rho_l, rho_g, mu_l, mu_g, sigma, d, g, grav=9.81):
    C0 = 1+0.12*(1-x)
    u_gj = 1.18*(1-x)*(grav*sigma*(rho_l-rho_g)/(rho_l**2))**0.25
    alpha = (x/rho_g)/(C0*((1-x)/rho_l+x/rho_g)+(u_gj/g))
    return alpha


def void_fraction_correlations(
    df: pd.DataFrame,
    *,
    ggrav: float = 9.78,
) -> pd.DataFrame:
    """Apply void-fraction correlations row-wise to a properties dataframe.

    Parameters
    ----------
    df
        Output of :func:`~src.combined_properties.properties_df`, with at least
        ``x``, ``rhol``, ``rhog``, ``mul``, ``mug``, ``sigma``, ``g``, ``d``.
    ggrav
        Gravitational acceleration [m/s²] passed to correlations that use it.

    Returns
    -------
    Copy of ``df`` with void-fraction correlation columns appended.
    """
    required = {"x", "rhol", "rhog", "mul", "mug", "sigma", "g", "d"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    out = df.copy()

    out["alpha_h"] = [
        alpha_homo(xi, rl, rg) for xi, rl, rg in zip(out["x"], out["rhol"], out["rhog"])
    ]
    out["lockhart"] = [
        alpha_lockhart_martinelli_1949(xi, rl, rg, mul, mug)
        for xi, rl, rg, mul, mug in zip(out["x"], out["rhol"], out["rhog"], out["mul"], out["mug"])
    ]
    out["zivi"] = [
        alpha_zivi_1964(xi, rl, rg) for xi, rl, rg in zip(out["x"], out["rhol"], out["rhog"])
    ]
    out["premoli"] = [
        alpha_premoli_1970(xi, rl, rg, mul, sigma, g, di)
        for xi, rl, rg, mul, sigma, g, di in zip(
            out["x"], out["rhol"], out["rhog"], out["mul"], out["sigma"], out["g"], out["d"]
        )
    ]
    out["kanizawa"] = [
        alpha_kanizawa_ribastki_2015(xi, rl, rg, mul, mug, g, di, ggrav)
        for xi, rl, rg, mul, mug, g, di in zip(
            out["x"], out["rhol"], out["rhog"], out["mul"], out["mug"], out["g"], out["d"]
        )
    ]
    out["tibirica"] = [
        alpha_tibirica_2017(xi, rl, rg, g, di, ggrav)
        for xi, rl, rg, g, di in zip(out["x"], out["rhol"], out["rhog"], out["g"], out["d"])
    ]
    out["hughmark"] = [
        alpha_hughmark_1962(xi, rl, rg, mul, mug, di, g, ggrav)
        for xi, rl, rg, mul, mug, di, g in zip(
            out["x"], out["rhol"], out["rhog"], out["mul"], out["mug"], out["d"], out["g"]
        )
    ]
    out["rouhani"] = [
        alpha_rouhani_1970(xi, rl, rg, mul, mug, sigma, di, g, ggrav)
        for xi, rl, rg, mul, mug, sigma, di, g in zip(
            out["x"],
            out["rhol"],
            out["rhog"],
            out["mul"],
            out["mug"],
            out["sigma"],
            out["d"],
            out["g"],
        )
    ]
    out["steiner"] = [
        alpha_steiner_1993(xi, rl, rg, mul, mug, sigma, di, g, ggrav)
        for xi, rl, rg, mul, mug, sigma, di, g in zip(
            out["x"],
            out["rhol"],
            out["rhog"],
            out["mul"],
            out["mug"],
            out["sigma"],
            out["d"],
            out["g"],
        )
    ]

    return out