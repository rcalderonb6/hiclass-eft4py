"""Tests for tabulated inputs (gravity_model = alphas_table, expansion_model = table),
the wrapper's array/length handling, and failed-shooting reporting.

Run from this directory, with hiclassy installed:

    python -m pytest test_tables_smg.py -q --confcutdir=.

(--confcutdir=. stops pytest from importing the repository root, whose
__init__.py belongs to the installed hiclassy package.)

The closure tests (V1) run a model covariantly, feed its own background and
alphas back through the table models, and compare the spectra.
"""

import os
import re
import sys

import numpy as np
import pytest

hiclassy = pytest.importorskip("hiclassy")
from hiclassy import Class, CosmoComputationError, CosmoSevereError, HiClass  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from tables_smg import background_tables, table_parameters, write_table_files  # noqa: E402

# Converged shooting: hi_class's default tolerances leave e.g. nKGB's budget
# unclosed at the 1e-3 level, which would dominate the comparison.
SHOOTING = {"tol_shooting_deltax_rel": 1e-12, "tol_shooting_deltaF": 1e-10}
SPECTRA = {"output": "tCl,pCl,lCl,mPk", "lensing": "yes", "l_max_scalars": 2600,
           "P_k_max_1/Mpc": 2.0, "z_pk": "0,1", "output_background_smg": 10, **SHOOTING}
MODEL_KEYS = {"gravity_model", "parameters_smg", "expansion_model", "expansion_smg",
              "M2_tuning_smg", "M2_today_smg", "tuning_index_smg", "tuning_dxdy_guess_smg",
              "param_shoot_M2_smg", "attractor_ic_smg", "Omega_smg"}
IGNORED = re.compile(r"^(root|overwrite_root|output|write .*|.*_verbose|output_background_smg|"
                     r"k_output_values|headers|format|write_.*)$")

# Targets (dev plan, V1)
TARGET_SIGMA_CV = 0.01
TARGET_PK = 1e-3
TARGET_SCALAR = 1e-6


def read_ini(path):
    pars = {}
    with open(path) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if "=" in line:
                key, value = [x.strip() for x in line.split("=", 1)]
                if not IGNORED.match(key):
                    pars[key] = value
    return pars


def model(name):
    if name == "propto_omega":
        return {**read_ini(os.path.join(ROOT, "hi_class.ini")), "parameters_smg": "1., 0.625, 0.3, 0., 1."}
    return read_ini(os.path.join(ROOT, "gravity_models", name + ".ini"))


def compute(pars):
    c = Class()
    c.set(pars)
    c.compute()
    return c


def observables(c):
    cl = c.lensed_cl(2500)
    k = np.logspace(-4, 0, 100)
    pk = {z: np.array([c.pk(kk, z) for kk in k]) for z in (0.0, 1.0)}
    return {"cl": cl, "pk": pk, "h": c.h(), "theta": c.theta_s_100()}


def worst_binned_sigma(ref, other, key, dl=30):
    """Worst |sum of differences| over the cosmic-variance error of the bin sum."""
    ell = np.arange(2501)
    a, b = ref["cl"][key], other["cl"][key]
    if key == "te":
        sigma = np.sqrt((ref["cl"]["tt"] * ref["cl"]["ee"] + a**2) / (2 * ell + 1.0))
    else:
        sigma = np.sqrt(2.0 / (2 * ell + 1.0)) * np.abs(a)
    worst = 0.0
    for lo in range(2, 2501, dl):
        s = slice(lo, min(lo + dl, 2501))
        worst = max(worst, abs(np.sum(b[s] - a[s])) / np.sqrt(np.sum(sigma[s] ** 2)))
    return worst


def table_inputs(pars, cov, **kwargs):
    base = {k: v for k, v in pars.items() if k not in MODEL_KEYS}
    tabs = background_tables(cov.get_background(), kineticity_safe=float(pars.get("kineticity_safe_smg", 0.0)), **kwargs)
    return base, tabs


@pytest.fixture(scope="module")
def propto_omega_tables():
    pars = {**model("propto_omega"), "output": "", **SHOOTING}
    cov = compute({**pars, "output_background_smg": 10})
    base, tabs = table_inputs(pars, cov, n_nodes=200)
    return base, tabs


# ---------------------------------------------------------------------------
# Wrapper
# ---------------------------------------------------------------------------

def test_class_alias():
    assert Class is HiClass


def test_array_values_round_trip_exactly():
    from hiclassy._hiclassy import _format_parameter_value
    x = np.random.default_rng(1).normal(size=3000) * 10.0 ** np.random.default_rng(2).integers(-300, 300, 3000)
    s = _format_parameter_value(x)
    assert np.array_equal(np.array([float(v) for v in s.split(",")]), x)
    assert _format_parameter_value([1.0, 2.5]) == "1.0,2.5"
    assert _format_parameter_value("a, b") == "a, b"
    assert _format_parameter_value(0.5) == 0.5


def test_oversized_value_raises():
    c = Class()
    c.set({"output": "", "some_long_list": np.full(20000, 0.1234567890123)})
    with pytest.raises(CosmoSevereError, match="characters long"):
        c.compute()


# ---------------------------------------------------------------------------
# Table input
# ---------------------------------------------------------------------------

def test_file_and_list_inputs_agree(propto_omega_tables, tmp_path):
    base, tabs = propto_omega_tables
    lists = compute({**base, **table_parameters(tabs), "output": "tCl"})
    write_table_files(tabs, tmp_path / "a.dat", tmp_path / "e.dat")
    files = compute({**base, "gravity_model": "alphas_table", "expansion_model": "table",
                     "Omega_smg": -1, "Omega_Lambda": 0, "Omega_fld": 0,
                     "alphas_table_file_smg": str(tmp_path / "a.dat"),
                     "expansion_table_file_smg": str(tmp_path / "e.dat"), "output": "tCl"})
    assert np.array_equal(lists.raw_cl(2000)["tt"], files.raw_cl(2000)["tt"])


def test_table_mode_closes_budget_without_shooting(propto_omega_tables):
    base, tabs = propto_omega_tables
    c = compute({**base, **table_parameters(tabs), "input_verbose": 1})
    bg = c.get_background()
    H0 = c.h() * 100 / 299792.458
    # exact up to the accuracy of CLASS's background integration (~1e-9)
    assert abs(bg["H [1/Mpc]"][-1] / H0 - 1) < 1e-8


def _broken(tabs, **changes):
    t = {k: np.array(v, copy=True) for k, v in tabs.items()}
    for k, f in changes.items():
        t[k] = f(t[k])
    return t


@pytest.mark.parametrize("changes, message", [
    ({"lna": lambda x: x[::-1]}, "strictly increasing"),
    ({"lna": lambda x: x - 0.1}, "last node"),
    ({"alpha_M": lambda x: 2.0 * x + 0.1}, "alpha_M is inconsistent"),
    ({"Delta_M2": lambda x: x - 2.0}, "M_\\*\\^2"),
    ({"rho_de": lambda x: 1.1 * x}, "normalised"),
    ({"alpha_K": lambda x: np.where(np.arange(len(x)) == 5, np.nan, x)}, "non-finite|invalid number"),
])
def test_invalid_tables_raise(propto_omega_tables, changes, message):
    base, tabs = propto_omega_tables
    c = Class()
    c.set({**base, **table_parameters(_broken(tabs, **changes))})
    with pytest.raises((CosmoSevereError, CosmoComputationError), match=message):
        c.compute()


def test_mismatched_list_lengths_raise(propto_omega_tables):
    base, tabs = propto_omega_tables
    params = table_parameters(tabs)
    params["table_alpha_B_smg"] = params["table_alpha_B_smg"][:-1]
    c = Class()
    c.set({**base, **params})
    with pytest.raises(CosmoSevereError, match="same length"):
        c.compute()


def test_file_and_lists_together_raise(propto_omega_tables, tmp_path):
    base, tabs = propto_omega_tables
    write_table_files(tabs, tmp_path / "a.dat", tmp_path / "e.dat")
    c = Class()
    c.set({**base, **table_parameters(tabs), "alphas_table_file_smg": str(tmp_path / "a.dat")})
    with pytest.raises(CosmoSevereError, match="not both"):
        c.compute()


def test_malformed_file_raises(propto_omega_tables, tmp_path):
    base, tabs = propto_omega_tables
    write_table_files(tabs, tmp_path / "a.dat", tmp_path / "e.dat")
    with open(tmp_path / "a.dat", "a") as f:
        f.write("0.5 1 2 x 4 5\n")
    params = {k: v for k, v in table_parameters(tabs).items() if not k.startswith("table_") or "_de_" in k}
    c = Class()
    c.set({**base, **params, "alphas_table_file_smg": str(tmp_path / "a.dat")})
    with pytest.raises(CosmoSevereError, match="cannot read a number"):
        c.compute()


# ---------------------------------------------------------------------------
# V1: closure inside hi_class
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", ["propto_omega", "quintessence_monomial", "nkgb", "brans_dicke"])
def test_v1_closure(name):
    pars = {**model(name), **SPECTRA}
    cov = compute(pars)
    base, tabs = table_inputs(pars, cov, n_nodes=600, a_min=1e-8)
    tab = compute({**base, **table_parameters(tabs)})
    ref, out = observables(cov), observables(tab)
    for key in ("tt", "te", "ee", "pp"):
        assert worst_binned_sigma(ref, out, key) < TARGET_SIGMA_CV, key
    for z in (0.0, 1.0):
        assert np.max(np.abs(out["pk"][z] / ref["pk"][z] - 1)) < TARGET_PK
    assert abs(out["h"] / ref["h"] - 1) < TARGET_SCALAR
    assert abs(out["theta"] / ref["theta"] - 1) < TARGET_SCALAR


# ---------------------------------------------------------------------------
# Tables that start late: early-time continuation
# ---------------------------------------------------------------------------

def _closure(name, a_min, **extra):
    pars = {**model(name), **SPECTRA}
    cov = compute(pars)
    base, tabs = table_inputs(pars, cov, n_nodes=600, a_min=a_min)
    tab = compute({**base, **table_parameters(tabs), **extra})
    ref, out = observables(cov), observables(tab)
    return (max(worst_binned_sigma(ref, out, k) for k in ("tt", "te", "ee", "pp")),
            max(np.max(np.abs(out["pk"][z] / ref["pk"][z] - 1)) for z in (0.0, 1.0)),
            abs(out["theta"] / ref["theta"] - 1))


@pytest.mark.parametrize("name, a_min", [
    ("brans_dicke", 1e-6),          # rho_smg tracks radiation; continued via rho_phi
    ("quintessence_monomial", 1e-5),
    ("nkgb", 1e-4),                 # alpha ~ rho_bg^(-4/3) recovered across equality
    ("propto_omega", 1e-4),         # alpha ~ Omega_smg ~ rho_bg^(-1)
])
def test_late_start_closure(name, a_min):
    cl, pk, theta = _closure(name, a_min)
    assert cl < TARGET_SIGMA_CV and pk < TARGET_PK and theta < TARGET_SCALAR


def test_undetermined_continuation_raises():
    pars = {**model("quintessence_monomial"), "output": "", **SHOOTING}
    cov = compute({**pars, "output_background_smg": 10})
    base, tabs = table_inputs(pars, cov, n_nodes=600, a_min=1e-4)
    c = Class()
    c.set({**base, **table_parameters(tabs)})
    with pytest.raises(CosmoComputationError, match="not determined by the table"):
        c.compute()


# ---------------------------------------------------------------------------
# Failed shootings are reported
# ---------------------------------------------------------------------------

def test_failed_shooting_raises():
    pars = {**model("brans_dicke"), "M2_today_smg": 0.01, "output": ""}
    c = Class()
    c.set(pars)
    with pytest.raises(CosmoComputationError, match="Shooting failed"):
        c.compute()
