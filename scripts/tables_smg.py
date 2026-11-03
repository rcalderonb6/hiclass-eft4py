"""Build ``alphas_table`` / ``expansion_model = table`` inputs for hi_class.

The tables are made from a hi_class background (``HiClass.get_background()``).
The closure test (V1) uses this: a model is run covariantly, then its own
background and alphas are fed back through the table models, and the spectra
are compared.

Columns are in hi_class's own conventions, i.e. those of the background output:
``kineticity_smg``, ``braiding_smg``, ``M2_running_smg``, ``tensor_excess_smg`` and
``M*^2_smg``. ``rho_smg`` is normalised to 1 at a = 1.

Example::

    from hiclassy import Class
    from tables_smg import table_parameters

    cov = Class(); cov.set(params); cov.compute()
    tab = Class(); tab.set({**base_params, **table_parameters(cov.get_background())})
"""

import numpy as np

TABLE_KEYS = (
    "table_lna_smg", "table_alpha_K_smg", "table_alpha_B_smg",
    "table_alpha_M_smg", "table_alpha_T_smg", "table_Delta_M2_smg",
)
EXPANSION_KEYS = ("table_lna_de_smg", "table_rho_de_smg", "table_w_de_smg")


def background_tables(background, n_nodes=600, a_min=1e-8, kineticity_safe=0.0):
    """Resample a hi_class background into table columns on a uniform ln a grid.

    Args:
        background: The dictionary returned by `HiClass.get_background()`.
        n_nodes: Number of nodes, uniform in ln a on [ln a_min, 0].
        a_min: Scale factor of the first node.
        kineticity_safe: The `kineticity_safe_smg` of the run that produced
            `background`. hi_class adds it to the kineticity it outputs, so it
            is subtracted here; set the same value in the table run.

    Returns:
        A dictionary of 1D arrays: `lna`, `alpha_K`, `alpha_B`, `alpha_M`,
        `alpha_T`, `Delta_M2`, `rho_de` (normalised to 1 at a = 1) and `w_de`.
    """
    from scipy.interpolate import CubicSpline

    lna_bg = -np.log1p(background["z"])
    order = np.argsort(lna_bg)
    lna_bg = lna_bg[order]
    keep = np.concatenate(([True], np.diff(lna_bg) > 0))
    lna_bg = lna_bg[keep]

    def column(name):
        return np.asarray(background[name])[order][keep]

    lna = np.linspace(np.log(a_min), 0.0, n_nodes)
    if lna[0] < lna_bg[0]:
        raise ValueError("a_min = %g is earlier than the background table" % a_min)

    rho = column("(.)rho_smg")
    p = column("(.)p_smg")
    log_rho = CubicSpline(lna_bg, np.log(rho))(lna)
    out = {
        "lna": lna,
        "alpha_K": CubicSpline(lna_bg, column("kineticity_smg") - kineticity_safe)(lna),
        "alpha_B": CubicSpline(lna_bg, column("braiding_smg"))(lna),
        "alpha_M": CubicSpline(lna_bg, column("M2_running_smg"))(lna),
        "alpha_T": CubicSpline(lna_bg, column("tensor_excess_smg"))(lna),
        "Delta_M2": CubicSpline(lna_bg, column("D_M*^2_smg"))(lna),
        "rho_de": np.exp(log_rho - log_rho[-1]),
        "w_de": CubicSpline(lna_bg, p / rho)(lna),
    }
    return out


def table_parameters(tables, with_w=True):
    """Turn table columns into hi_class input parameters.

    Args:
        tables: Columns as returned by `background_tables`.
        with_w: Whether to pass `w` (then p = w rho exactly) or let hi_class
            derive it from d ln rho / d ln a.

    Returns:
        A parameter dictionary for `HiClass.set`, with `gravity_model =
        alphas_table`, `expansion_model = table` and `Omega_smg = -1`.
    """
    params = {
        "gravity_model": "alphas_table",
        "expansion_model": "table",
        "Omega_smg": -1,
        "Omega_Lambda": 0,
        "Omega_fld": 0,
    }
    for key, name in zip(TABLE_KEYS, ("lna", "alpha_K", "alpha_B", "alpha_M", "alpha_T", "Delta_M2")):
        params[key] = np.asarray(tables[name])
    params["table_lna_de_smg"] = np.asarray(tables["lna"])
    params["table_rho_de_smg"] = np.asarray(tables["rho_de"])
    if with_w:
        params["table_w_de_smg"] = np.asarray(tables["w_de"])
    return params


def write_table_files(tables, alphas_path, expansion_path, with_w=True):
    """Write the tables as files for `alphas_table_file_smg` / `expansion_table_file_smg`."""
    np.savetxt(alphas_path,
               np.column_stack([tables[k] for k in ("lna", "alpha_K", "alpha_B", "alpha_M", "alpha_T", "Delta_M2")]),
               header="ln a, alpha_K, alpha_B, alpha_M, alpha_T, Delta M_*^2 (hi_class conventions)",
               fmt="%.17g")
    cols = ["lna", "rho_de"] + (["w_de"] if with_w else [])
    np.savetxt(expansion_path, np.column_stack([tables[k] for k in cols]),
               header="ln a, rho_smg/rho_smg(a=1)" + (", w_smg" if with_w else ""), fmt="%.17g")
