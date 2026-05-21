import pandas as pd
from CoolProp.CoolProp import PropsSI

from .alpha_correlations import alpha_homo
from .dimensionless_number import Bond, Reynolds


def rho_x(x, rho_l, rho_g):
    rho = ((1 - x) / rho_l + x / rho_g) ** -1
    return rho


def rho_alpha(x, rho_l, rho_g):
    alpha_h = alpha_homo(x, rho_l, rho_g)
    rho = ((1 - alpha_h) / rho_l + alpha_h / rho_g) ** -1
    return rho


def Jv(g, x, d, rho_l, rho_g, ggrav=9.78):
    J = g * x / ((ggrav * d * rho_g * (rho_l - rho_g)) ** 0.5)
    return J


def Gw(rho_l, rho_g, d, sigma, ggrav=9.78):
    bd = Bond(d, rho_l, rho_g, sigma)
    G = rho_l * (ggrav * d) ** 0.5 * (0.54 - 0.96 * bd**-2 - 4.2 * bd**-1)
    return G


def properties_df(
    df: pd.DataFrame,
    *,
    d: float = 8e-3,
    fluid_l: str = "water",
    fluid_g: str = "air",
) -> pd.DataFrame:
    """Compute test-section thermophysical properties and flow parameters.

    Parameters
    ----------
    df
        Experimental runs with at least ``run_id``, ``Ttp`` [°C], ``Ptp`` [Pa],
        ``jl`` and ``jg`` [m/s] (superficial volumetric velocities).
    d
        Tube inner diameter [m].
    fluid_l, fluid_g
        CoolProp fluid names for liquid and gas phases.

    Returns
    -------
    DataFrame with ``run_id``, flow conditions, thermophysical properties at
    (``Ptp``, ``Ttp``), mass flux ``g``, quality ``x``, and phase Reynolds numbers.
    """
    required = {"run_id", "Ttp", "Ptp", "jl", "jg"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    out = df.loc[:, ["run_id", "Ttp", "Ptp", "jl", "jg"]].copy()
    out["d"] = d

    p = out["Ptp"].to_numpy()
    t = (out["Ttp"] + 273.15).to_numpy()

    out["mul"] = [PropsSI("V", "P", pi, "T", ti, fluid_l) for pi, ti in zip(p, t)]
    out["mug"] = [PropsSI("V", "P", pi, "T", ti, fluid_g) for pi, ti in zip(p, t)]
    out["rhol"] = [PropsSI("D", "P", pi, "T", ti, fluid_l) for pi, ti in zip(p, t)]
    out["rhog"] = [PropsSI("D", "P", pi, "T", ti, fluid_g) for pi, ti in zip(p, t)]
    out["sigma"] = [PropsSI("I", "T", ti, "Q", 0, fluid_l) for ti in t]

    out["Gl"] = out["jl"] * out["rhol"]
    out["Gg"] = out["jg"] * out["rhog"]
    out["g"] = out["Gl"] + out["Gg"]
    out["x"] = out["Gg"] / out["g"]

    out["Rel"] = Reynolds(out["g"] * (1 - out["x"]), out["d"], out["mul"])
    out["Reg"] = Reynolds(out["g"] * out["x"], out["d"], out["mug"])

    return out
