#!/usr/bin/env python3
"""
PFAS Censoring-Aware Screening Tool (PFAS-CAST), Python screening engine v1.1.1
==============================================================================

Companion to:
  Kharel, G., Kirilenko, A. P., Harvey, O. R. and Zou, C. B. Censoring-aware models quantify systematic underestimation of PFAS occurrence in US drinking-water monitoring. Cell Reports Sustainability (2026). https://doi.org/10.1016/j.crsus.2026.100893

Local command-line companion to the browser-based PFAS Censoring-Aware Screening
Tool (PFAS-CAST). This script implements the manuscript-aligned primary workflow used for
screening uploaded PFAS monitoring data against configurable benchmark
thresholds:

  • record-level left-censoring;
  • analyte-specific censored-lognormal likelihoods;
  • tier-aware parameter uncertainty (A, B2, B1, C);
  • record-level multiple imputation followed by repeated-record aggregation;
  • benchmark-ratio screening under zero substitution and multiple imputation;
  • Excel-compatible output for transparent review/audit.

Important scope statement
-------------------------
This is a screening and prioritization tool, not a regulatory-compliance tool.
It estimates benchmark-ratio shifts for uploaded monitoring records; it does not
reproduce the manuscript's national monitored-CWS population-weighted occurrence
estimand, delivered-load metrics, intake dose, or toxicity-weighted risk.
Default benchmarks reflect the April 2024 PFAS NPDWR framework as configurable
screening thresholds. Users should verify current regulatory status and local
requirements before applying any benchmark for decision support.

The tool intentionally implements the manuscript's primary marginal MI engine.
It does not include the manuscript's secondary copula diagnostic, which is not
load-bearing in the paper.

Input files may be comma- or tab-delimited (EPA UCMR 5 Data Finder exports are tab-delimited).
UCMR 5 column names, including AnalyticalResultValue(µg/L), AnalyticalResultsSign and
SamplePointID, are recognized automatically.

Dependencies: numpy, pandas, scipy, openpyxl
License: MIT
"""

from __future__ import annotations

import argparse
import dataclasses as dc
import hashlib
import json
import math
import os
import re
import sys
import time
import warnings
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import norm

try:
    import openpyxl  # noqa: F401  # required by pandas ExcelWriter
except Exception as exc:  # pragma: no cover
    raise RuntimeError("openpyxl is required for .xlsx output. Install with: pip install openpyxl") from exc

TOOL_VERSION = "1.1.1"
INTERNAL_UNIT = "ug/L"
EPS = 1e-15

DEFAULT_INDIVIDUAL_BENCHMARKS_UG_L = {
    "PFOA": 0.004,
    "PFOS": 0.004,
    "PFHXS": 0.010,
    "PFNA": 0.010,
    "HFPO-DA": 0.010,
}

DEFAULT_HI_DENOMINATORS_UG_L = {
    "PFHXS": 0.010,
    "PFNA": 0.010,
    "HFPO-DA": 0.010,
    "PFBS": 2.000,
}

COLUMN_ALIASES = {
    "system_id": [
        "pwsid", "pws_id", "pws", "system_id", "systemid", "water_system_id",
        "publicwatersystemid", "public_water_system_id", "pws identification number",
    ],
    "location_id": [
        "samplinglocationid", "sampling_location_id", "location_id", "locationid",
        "samplepointid", "sample_point_id", "sample_point", "sampling_point",
        "facility_id", "entry_point", "entrypoint",
        "samplelocation", "samplinglocation",
    ],
    "analyte": [
        "analyte", "contaminant", "contaminantname", "contaminant_name",
        "parameter", "chemical", "compound", "component",
    ],
    "value": [
        "analyticalresultvalue", "analyticalresultvaluegl", "analyticalresultvalueugl",
        "analytical_result_value", "result", "result_value",
        "concentration", "value", "detect_value", "measured_value", "reported_value",
    ],
    "detect": [
        "detect", "detected", "detect_flag", "detection_flag", "is_detect",
        "analyticalresultssign", "analyticalresultsign", "result_sign", "sign",
        "result_status", "resultqualifier", "result_qualifier", "detectindicator",
    ],
    "reporting_limit": [
        "mrl", "minimumreportinglevel", "minimum_reporting_level", "reporting_limit",
        "rl", "mdl", "lod", "limit", "censoring_limit", "detection_limit",
    ],
}

ANALYTE_ALIASES = {
    "PFOA": ["PFOA", "PERFLUOROOCTANOIC ACID"],
    "PFOS": ["PFOS", "PERFLUOROOCTANESULFONIC ACID", "PERFLUOROOCTANE SULFONIC ACID"],
    "PFHXS": ["PFHXS", "PFHXS ACID", "PERFLUOROHEXANESULFONIC ACID", "PERFLUOROHEXANE SULFONIC ACID"],
    "PFNA": ["PFNA", "PERFLUORONONANOIC ACID"],
    "HFPO-DA": ["HFPO-DA", "HFPO DA", "HFPO_DA", "GENX", "GENX CHEMICAL", "HFPO DIMER ACID"],
    "PFBS": ["PFBS", "PERFLUOROBUTANESULFONIC ACID", "PERFLUOROBUTANE SULFONIC ACID"],
    "PFBA": ["PFBA", "PERFLUOROBUTANOIC ACID"],
    "PFHPA": ["PFHPA", "PERFLUOROHEPTANOIC ACID"],
    "PFHXA": ["PFHXA", "PERFLUOROHEXANOIC ACID"],
    "PFDA": ["PFDA", "PERFLUORODECANOIC ACID"],
    "PFUNA": ["PFUNA", "PFUNDA", "PERFLUOROUNDECANOIC ACID"],
}


def _norm_key(x: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())


def _norm_analyte(x: Any) -> str:
    s = str(x).strip().upper()
    s2 = re.sub(r"[^A-Z0-9]+", "", s)
    for canonical, aliases in ANALYTE_ALIASES.items():
        for alias in aliases:
            if s2 == re.sub(r"[^A-Z0-9]+", "", alias.upper()):
                return canonical
    return s.replace("_", "-")


def _parse_number(x: Any) -> float:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return np.nan
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x)
    s = str(x).strip().replace(",", "")
    if not s:
        return np.nan
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    return float(m.group(0)) if m else np.nan


def _parse_detect_flag(flag: Any, raw_value: Any = None) -> Optional[bool]:
    if flag is None or (isinstance(flag, float) and np.isnan(flag)):
        if isinstance(raw_value, str) and raw_value.strip().startswith("<"):
            return False
        return None
    if isinstance(flag, (bool, np.bool_)):
        return bool(flag)
    if isinstance(flag, (int, float, np.integer, np.floating)) and not np.isnan(flag):
        return bool(int(flag))
    s = str(flag).strip().lower()
    if not s:
        if isinstance(raw_value, str) and raw_value.strip().startswith("<"):
            return False
        return None
    true_set = {"1", "y", "yes", "true", "t", "detected", "detect", "d", "present", "positive", "hit", "=", ">"}
    false_set = {"0", "n", "no", "false", "f", "nondetect", "non-detect", "non detect", "nd", "<mrl", "below", "censored"}
    if s in true_set:
        return True
    if s in false_set or s.startswith("<"):
        return False
    if isinstance(raw_value, str) and raw_value.strip().startswith("<"):
        return False
    return None


def _unit_factor_to_ug_l(unit: str) -> float:
    u = unit.strip().lower().replace("μ", "µ")
    if u in {"ug/l", "µg/l", "microg/l", "microgram/l", "micrograms/l"}:
        return 1.0
    if u in {"ng/l", "ppt", "nanog/l", "nanogram/l", "nanograms/l"}:
        return 0.001
    raise ValueError(f"Unsupported unit '{unit}'. Use ug/L or ng/L/ppt.")


def _find_column(df: pd.DataFrame, key: str, explicit: Optional[str] = None) -> Optional[str]:
    if explicit:
        if explicit not in df.columns:
            raise ValueError(f"Explicit column '{explicit}' not found for {key}.")
        return explicit
    norm_to_col = {_norm_key(c): c for c in df.columns}
    for alias in COLUMN_ALIASES[key]:
        k = _norm_key(alias)
        if k in norm_to_col:
            return norm_to_col[k]
    return None


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


@dc.dataclass
class FitResult:
    analyte: str
    n_total: int
    n_detect: int
    n_censored: int
    tier: str
    mode: str
    mu: float
    sigma: float
    eta: float
    loglik: float
    success: bool
    warnings: List[str] = dc.field(default_factory=list)
    cov: Optional[np.ndarray] = None
    draw_bank: Optional[np.ndarray] = None
    grid_params: Optional[np.ndarray] = None
    grid_weights: Optional[np.ndarray] = None
    eb_mu_bar: Optional[float] = None
    eb_s2: Optional[float] = None
    eb_sigma_bar: Optional[float] = None

    def draw_params(self, rng: np.random.Generator) -> Tuple[float, float]:
        if self.tier == "A" and self.cov is not None and np.all(np.isfinite(self.cov)):
            for _ in range(50):
                theta = rng.multivariate_normal(np.array([self.mu, self.eta]), self.cov)
                sigma = float(np.exp(theta[1]))
                if np.isfinite(theta[0]) and 0.02 <= sigma <= 10:
                    return float(theta[0]), sigma
            return self.mu, self.sigma
        if self.tier == "B2" and self.draw_bank is not None and len(self.draw_bank) > 0:
            theta = self.draw_bank[int(rng.integers(0, len(self.draw_bank)))]
            return float(theta[0]), float(np.exp(theta[1]))
        if self.tier == "B1" and self.grid_params is not None and self.grid_weights is not None:
            idx = int(rng.choice(len(self.grid_weights), p=self.grid_weights))
            theta = self.grid_params[idx]
            return float(theta[0]), float(np.exp(theta[1]))
        return self.mu, self.sigma


def lognormal_loglik(mu: float, eta: float, x_det: np.ndarray, L_cens: np.ndarray) -> float:
    sigma = float(np.exp(eta))
    if not np.isfinite(mu) or not np.isfinite(sigma) or sigma <= 0:
        return -np.inf
    ll = 0.0
    if x_det.size:
        lx = np.log(x_det)
        ll += float(np.sum(norm.logpdf(lx, loc=mu, scale=sigma) - lx))
    if L_cens.size:
        z = (np.log(L_cens) - mu) / sigma
        lcdf = norm.logcdf(z)
        if not np.all(np.isfinite(lcdf)):
            lcdf = np.where(np.isfinite(lcdf), lcdf, -1e9)
        ll += float(np.sum(lcdf))
    return ll


def mle_point(x_det: np.ndarray, L_cens: np.ndarray, init: Optional[Tuple[float, float]] = None) -> Tuple[float, float, float, bool, str]:
    x_det = np.asarray(x_det, dtype=float)
    L_cens = np.asarray(L_cens, dtype=float)
    x_det = x_det[np.isfinite(x_det) & (x_det > 0)]
    L_cens = L_cens[np.isfinite(L_cens) & (L_cens > 0)]
    if init is None:
        if x_det.size >= 2:
            logs = np.log(x_det)
            mu0 = float(np.mean(logs))
            sigma0 = float(max(np.std(logs, ddof=1), 0.35))
        elif x_det.size == 1:
            mu0 = float(np.log(x_det[0]))
            sigma0 = 1.0
        elif L_cens.size:
            mu0 = float(np.log(np.median(L_cens) * 0.20))
            sigma0 = 1.0
        else:
            mu0, sigma0 = -7.0, 1.0
    else:
        mu0, sigma0 = init
    eta0 = float(np.log(max(sigma0, 0.05)))
    mu_bounds = (mu0 - 12, mu0 + 12)
    bounds = [mu_bounds, (np.log(0.03), np.log(8.0))]

    def obj(theta: np.ndarray) -> float:
        ll = lognormal_loglik(float(theta[0]), float(theta[1]), x_det, L_cens)
        if not np.isfinite(ll):
            return 1e100
        return -ll

    best = None
    for start in [np.array([mu0, eta0]), np.array([mu0 - 1.0, eta0]), np.array([mu0 + 1.0, eta0])]:
        res = minimize(obj, start, method="L-BFGS-B", bounds=bounds, options={"maxiter": 1000, "ftol": 1e-9})
        if best is None or res.fun < best.fun:
            best = res
    if best is None:
        return mu0, eta0, -np.inf, False, "optimizer_not_run"
    mu, eta = float(best.x[0]), float(best.x[1])
    return mu, eta, -float(best.fun), bool(best.success), str(best.message)


def numerical_hessian_nll(mu: float, eta: float, x_det: np.ndarray, L_cens: np.ndarray) -> Optional[np.ndarray]:
    theta = np.array([mu, eta], dtype=float)
    h = np.array([max(1e-4, abs(mu) * 1e-5), 1e-4], dtype=float)

    def nll(th: np.ndarray) -> float:
        ll = lognormal_loglik(float(th[0]), float(th[1]), x_det, L_cens)
        return -ll if np.isfinite(ll) else 1e100

    H = np.zeros((2, 2), dtype=float)
    f0 = nll(theta)
    if not np.isfinite(f0):
        return None
    for i in range(2):
        ei = np.zeros(2); ei[i] = h[i]
        H[i, i] = (nll(theta + ei) - 2 * f0 + nll(theta - ei)) / (h[i] ** 2)
    for i in range(2):
        for j in range(i + 1, 2):
            ei = np.zeros(2); ej = np.zeros(2)
            ei[i] = h[i]; ej[j] = h[j]
            H[i, j] = H[j, i] = (
                nll(theta + ei + ej) - nll(theta + ei - ej) - nll(theta - ei + ej) + nll(theta - ei - ej)
            ) / (4 * h[i] * h[j])
    if not np.all(np.isfinite(H)):
        return None
    try:
        vals, vecs = np.linalg.eigh(H)
        vals = np.maximum(vals, 1e-6)
        Hpd = (vecs * vals) @ vecs.T
        cov = np.linalg.inv(Hpd)
        vals2, vecs2 = np.linalg.eigh(cov)
        vals2 = np.clip(vals2, 1e-8, 4.0)
        return (vecs2 * vals2) @ vecs2.T
    except Exception:
        return None


def build_b2_bank(
    mu: float,
    eta: float,
    x_det: np.ndarray,
    L_det: np.ndarray,
    L_cens: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = 200,
) -> np.ndarray:
    sigma = float(np.exp(eta))
    bank: List[Tuple[float, float]] = []
    L_det = np.asarray(L_det, dtype=float)
    if L_det.size != x_det.size:
        L_det = np.zeros_like(x_det)
    for _ in range(n_boot):
        if x_det.size:
            c0 = np.zeros_like(x_det, dtype=float)
            ok = np.isfinite(L_det) & (L_det > 0)
            c0[ok] = norm.cdf((np.log(L_det[ok]) - mu) / sigma)
            c0 = np.clip(c0, 0.0, 1.0 - 1e-12)
            u = c0 + (1 - c0) * rng.random(x_det.size)
            u = np.clip(u, 1e-12, 1 - 1e-12)
            x_sim = np.exp(mu + sigma * norm.ppf(u))
        else:
            x_sim = x_det
        mu_b, eta_b, _, ok_b, _ = mle_point(x_sim, L_cens, init=(mu, sigma))
        if ok_b and np.isfinite(mu_b) and np.isfinite(eta_b):
            bank.append((mu_b, eta_b))
    if not bank:
        bank = [(mu, eta)]
    return np.asarray(bank, dtype=float)


def build_b1_grid(mu: float, eta: float, x_det: np.ndarray, L_cens: np.ndarray, ng: int = 40) -> Tuple[np.ndarray, np.ndarray]:
    sigma = float(np.exp(eta))
    mu_span = max(2.0, 2.5 * sigma)
    eta_span = 1.25
    mu_grid = np.linspace(mu - mu_span, mu + mu_span, ng)
    eta_grid = np.linspace(max(np.log(0.05), eta - eta_span), min(np.log(6.0), eta + eta_span), ng)
    params = []
    lls = []
    for m in mu_grid:
        for e in eta_grid:
            params.append((float(m), float(e)))
            lls.append(lognormal_loglik(float(m), float(e), x_det, L_cens))
    params_arr = np.asarray(params, dtype=float)
    lls_arr = np.asarray(lls, dtype=float)
    finite = np.isfinite(lls_arr)
    if not np.any(finite):
        return np.asarray([[mu, eta]], dtype=float), np.asarray([1.0], dtype=float)
    lls_arr = np.where(finite, lls_arr, np.max(lls_arr[finite]) - 1e6)
    w = np.exp(lls_arr - np.max(lls_arr))
    if not np.isfinite(w).all() or w.sum() <= 0:
        w = np.ones(len(w), dtype=float)
    w = w / w.sum()
    return params_arr, w


def fit_all_analytes(df: pd.DataFrame, seed: int, b2_boot: int = 200, b1_grid: int = 40) -> Dict[str, FitResult]:
    rng = np.random.default_rng(seed + 1009)
    fits: Dict[str, FitResult] = {}
    preliminary: Dict[str, Tuple[float, float, float, bool, str]] = {}

    for analyte, g in df.groupby("analyte"):
        det = g["detect"].to_numpy(dtype=bool)
        x_det = g.loc[det, "value_ug_l"].to_numpy(dtype=float)
        L_cens = g.loc[~det, "rl_ug_l"].to_numpy(dtype=float)
        x_det = x_det[np.isfinite(x_det) & (x_det > 0)]
        L_cens = L_cens[np.isfinite(L_cens) & (L_cens > 0)]
        if x_det.size > 0:
            preliminary[analyte] = mle_point(x_det, L_cens)

    tierA_mus, tierA_sigmas = [], []
    for analyte, (mu, eta, _, ok, _) in preliminary.items():
        n_det = int(df[(df["analyte"] == analyte) & (df["detect"])].shape[0])
        if ok and n_det >= 50:
            tierA_mus.append(mu)
            tierA_sigmas.append(float(np.exp(eta)))
    if tierA_mus:
        eb_mu_bar = float(np.mean(tierA_mus))
        eb_s2 = float(max(np.var(tierA_mus, ddof=1) if len(tierA_mus) > 1 else 1.0, 0.25))
        eb_sigma_bar = float(np.median(tierA_sigmas))
    else:
        all_mus = [v[0] for v in preliminary.values() if v[3]]
        all_sig = [float(np.exp(v[1])) for v in preliminary.values() if v[3]]
        eb_mu_bar = float(np.mean(all_mus)) if all_mus else -7.0
        eb_s2 = float(max(np.var(all_mus, ddof=1) if len(all_mus) > 1 else 1.0, 0.25)) if all_mus else 1.0
        eb_sigma_bar = float(np.median(all_sig)) if all_sig else 1.0

    for analyte, g in df.groupby("analyte"):
        det = g["detect"].to_numpy(dtype=bool)
        x_det = g.loc[det, "value_ug_l"].to_numpy(dtype=float)
        L_det = g.loc[det, "rl_ug_l"].to_numpy(dtype=float)
        L_cens = g.loc[~det, "rl_ug_l"].to_numpy(dtype=float)
        x_det = x_det[np.isfinite(x_det) & (x_det > 0)]
        L_cens = L_cens[np.isfinite(L_cens) & (L_cens > 0)]
        n_total = int(len(g)); n_det = int(x_det.size); n_cens = int(L_cens.size)
        warnings_list: List[str] = []
        if n_det >= 50:
            tier = "A"
        elif n_det >= 10:
            tier = "B2"
        elif n_det >= 1:
            tier = "B1"
        else:
            tier = "C"

        if tier == "C":
            sigma = eb_sigma_bar
            eta = float(np.log(sigma))
            if L_cens.size:
                def neg_pen(mu_arr: Any) -> float:
                    mu_scalar = float(np.atleast_1d(mu_arr)[0])
                    ll = lognormal_loglik(mu_scalar, eta, np.array([], dtype=float), L_cens)
                    pen = -0.5 * ((mu_scalar - eb_mu_bar) ** 2) / eb_s2
                    return -(ll + pen)
                res = minimize_scalar(lambda z: neg_pen([z]), bounds=(eb_mu_bar - 8, eb_mu_bar + 4), method="bounded")
                mu = float(res.x) if res.success else float(np.log(np.median(L_cens) * 0.05))
                ok = bool(res.success)
                msg = str(res.message)
            else:
                mu = eb_mu_bar
                ok = False
                msg = "no_detections_no_censored_limits"
                warnings_list.append("Tier C fallback had no censoring limits.")
            ll = lognormal_loglik(mu, eta, np.array([], dtype=float), L_cens)
            warnings_list.append("Tier C zero-detection analyte: empirical-Bayes regularized fixed point; no parameter draws.")
            fits[analyte] = FitResult(
                analyte, n_total, n_det, n_cens, tier, "tier_c_empirical_bayes",
                mu, sigma, eta, ll, ok, warnings_list, None, None, None, None,
                eb_mu_bar, eb_s2, eb_sigma_bar,
            )
            continue

        mu, eta, ll, ok, msg = mle_point(x_det, L_cens)
        sigma = float(np.exp(eta))
        if not ok:
            warnings_list.append(f"MLE optimizer warning: {msg}")
        cov = None; bank = None; grid = None; weights = None
        mode = "censored_mle"
        if tier == "A":
            cov = numerical_hessian_nll(mu, eta, x_det, L_cens)
            if cov is None:
                warnings_list.append("Hessian covariance unavailable; parameter draws fall back to point estimate.")
                mode = "censored_mle_hessian_fallback"
            else:
                mode = "censored_mle_hessian"
        elif tier == "B2":
            bank = build_b2_bank(mu, eta, x_det, L_det, L_cens, rng, n_boot=b2_boot)
            mode = "censored_mle_conditional_bootstrap"
        elif tier == "B1":
            grid, weights = build_b1_grid(mu, eta, x_det, L_cens, ng=b1_grid)
            mode = "censored_mle_likelihood_grid"
        fits[analyte] = FitResult(
            analyte, n_total, n_det, n_cens, tier, mode, mu, sigma, eta, ll, ok,
            warnings_list, cov, bank, grid, weights, eb_mu_bar, eb_s2, eb_sigma_bar,
        )
    return fits


def parse_input_csv(
    input_path: Path,
    unit: str,
    explicit_cols: Dict[str, Optional[str]],
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    with input_path.open("rb") as fh:
        head = fh.read(4096).decode("utf-8", errors="replace")
    sep = "\t" if head.count("\t") > head.count(",") else ","
    try:
        raw = pd.read_csv(input_path, sep=sep, encoding="utf-8", low_memory=False)
    except UnicodeDecodeError:
        raw = pd.read_csv(input_path, sep=sep, encoding="latin-1", low_memory=False)
    factor = _unit_factor_to_ug_l(unit)
    cols = {k: _find_column(raw, k, explicit_cols.get(k)) for k in COLUMN_ALIASES}
    required = ["system_id", "analyte", "value"]
    missing = [k for k in required if cols[k] is None]
    if missing:
        raise ValueError(f"Required column(s) not found: {missing}. Use explicit --*-col arguments.")
    if cols["location_id"] is None:
        raw["__location_id"] = raw[cols["system_id"]].astype(str)
        cols["location_id"] = "__location_id"
    if cols["reporting_limit"] is None:
        raw["__reporting_limit"] = np.nan
        cols["reporting_limit"] = "__reporting_limit"
    if cols["detect"] is None:
        raw["__detect"] = None
        cols["detect"] = "__detect"

    records = []
    qa = defaultdict(int)
    unrec_detect = 0
    for idx, r in raw.iterrows():
        system = str(r[cols["system_id"]]).strip()
        loc = str(r[cols["location_id"]]).strip()
        analyte = _norm_analyte(r[cols["analyte"]])
        raw_val = r[cols["value"]]
        val = _parse_number(raw_val)
        rl = _parse_number(r[cols["reporting_limit"]])
        detect = _parse_detect_flag(r[cols["detect"]], raw_val)
        if detect is None:
            detect = not (isinstance(raw_val, str) and raw_val.strip().startswith("<"))
            if r[cols["detect"]] is not None and not (isinstance(r[cols["detect"]], float) and np.isnan(r[cols["detect"]])):
                unrec_detect += 1
        if isinstance(raw_val, str) and raw_val.strip().startswith("<") and (not np.isfinite(rl) or rl <= 0) and np.isfinite(val):
            rl = val
        if not system or system.lower() in {"nan", "none"} or not analyte or analyte.lower() == "nan":
            qa["dropped_missing_system_or_analyte"] += 1
            continue
        if detect:
            if not np.isfinite(val) or val <= 0:
                qa["dropped_detect_missing_positive_value"] += 1
                continue
        else:
            if not np.isfinite(rl) or rl <= 0:
                qa["dropped_censored_missing_positive_rl"] += 1
                continue
        records.append({
            "row_id": int(idx),
            "system_id": system,
            "location_id": loc if loc and loc.lower() not in {"nan", "none"} else system,
            "analyte": analyte,
            "detect": bool(detect),
            "value_ug_l": float(val * factor) if detect else np.nan,
            "rl_ug_l": float(rl * factor) if np.isfinite(rl) and rl > 0 else np.nan,
            "raw_value": raw_val,
        })
    df = pd.DataFrame(records)
    if df.empty:
        raise ValueError("No valid records after parsing input CSV.")
    qa.update({
        "input_rows": int(len(raw)),
        "valid_rows": int(len(df)),
        "dropped_rows": int(len(raw) - len(df)),
        "unrecognized_detect_flags": int(unrec_detect),
        "systems": int(df["system_id"].nunique()),
        "locations": int(df[["system_id", "location_id"]].drop_duplicates().shape[0]),
        "analytes": int(df["analyte"].nunique()),
        "detections": int(df["detect"].sum()),
        "censored": int((~df["detect"]).sum()),
        "unit_input": unit,
        "unit_internal": INTERNAL_UNIT,
        "columns": json.dumps(cols),
    })
    qa_df = pd.DataFrame([{"metric": k, "value": v} for k, v in qa.items()])
    return df, qa_df, dict(qa)


def load_benchmarks(path: Optional[Path]) -> Dict[str, Dict[str, float]]:
    if path is None:
        return {
            "individual_mcls_ug_l": DEFAULT_INDIVIDUAL_BENCHMARKS_UG_L.copy(),
            "hi_denominators_ug_l": DEFAULT_HI_DENOMINATORS_UG_L.copy(),
        }
    with path.open("r", encoding="utf-8") as f:
        obj = json.load(f)
    ind = obj.get("individual_mcls_ug_l") or obj.get("individual") or {}
    hi = obj.get("hi_denominators_ug_l") or obj.get("hazard_index") or {}
    return {
        "individual_mcls_ug_l": {_norm_analyte(k): float(v) for k, v in ind.items()},
        "hi_denominators_ug_l": {_norm_analyte(k): float(v) for k, v in hi.items()},
    }


def endpoint_ratios_from_location_means(location_means: pd.DataFrame, benchmarks: Dict[str, Dict[str, float]]) -> pd.DataFrame:
    rows = []
    analyte_cols = set(location_means.columns)
    for analyte, threshold in benchmarks["individual_mcls_ug_l"].items():
        if analyte in analyte_cols and threshold > 0:
            tmp = location_means[[analyte]].copy()
            tmp["endpoint"] = analyte
            tmp["ratio"] = tmp[analyte] / threshold
            rows.append(tmp[["endpoint", "ratio"]])
    hi_den = benchmarks.get("hi_denominators_ug_l", {})
    hi_terms = []
    for analyte, denom in hi_den.items():
        if analyte in analyte_cols and denom > 0:
            hi_terms.append(location_means[analyte] / denom)
    if hi_terms:
        hi = sum(hi_terms)
        rows.append(pd.DataFrame({"endpoint": "Hazard Index", "ratio": hi}, index=location_means.index))
    if not rows:
        return pd.DataFrame(columns=["system_id", "location_id", "endpoint", "ratio"])
    out = pd.concat(rows)
    out = out.reset_index()
    return out[["system_id", "location_id", "endpoint", "ratio"]]


def aggregate_records(df: pd.DataFrame, conc: np.ndarray) -> pd.DataFrame:
    tmp = df[["system_id", "location_id", "analyte"]].copy()
    tmp["conc"] = conc
    loc_mean = tmp.groupby(["system_id", "location_id", "analyte"], observed=True)["conc"].mean().unstack("analyte").fillna(0.0)
    return loc_mean


def impute_censored_values(df: pd.DataFrame, fits: Dict[str, FitResult], rng: np.random.Generator) -> np.ndarray:
    conc = np.where(df["detect"].to_numpy(dtype=bool), df["value_ug_l"].to_numpy(dtype=float), np.nan)
    analytes = df["analyte"].to_numpy()
    detect = df["detect"].to_numpy(dtype=bool)
    rl = df["rl_ug_l"].to_numpy(dtype=float)
    for analyte, fit in fits.items():
        idx = np.where((analytes == analyte) & (~detect))[0]
        if idx.size == 0:
            continue
        mu, sigma = fit.draw_params(rng)
        L = rl[idx]
        ok = np.isfinite(L) & (L > 0)
        draws = np.zeros(idx.size, dtype=float)
        if np.any(ok):
            cdf = norm.cdf((np.log(L[ok]) - mu) / sigma)
            cdf = np.clip(cdf, 1e-12, 1 - 1e-12)
            u = rng.random(np.sum(ok)) * cdf
            u = np.clip(u, 1e-12, 1 - 1e-12)
            draws[ok] = np.exp(mu + sigma * norm.ppf(u))
        if np.any(~ok):
            draws[~ok] = 0.0
        conc[idx] = draws
    return conc


def run_screening(
    df: pd.DataFrame,
    fits: Dict[str, FitResult],
    benchmarks: Dict[str, Dict[str, float]],
    m: int,
    seed: int,
    flag_rule: str = "mean",
    prob_threshold: float = 0.5,
    checkpoints: Optional[Sequence[int]] = None,
) -> Dict[str, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    zero_conc = np.where(df["detect"].to_numpy(dtype=bool), df["value_ug_l"].to_numpy(dtype=float), 0.0)
    zero_loc = aggregate_records(df, zero_conc)
    zero_endpoint = endpoint_ratios_from_location_means(zero_loc, benchmarks)
    if zero_endpoint.empty:
        raise ValueError("No uploaded analytes matched the benchmark endpoints.")
    zero_sys = zero_endpoint.groupby("system_id", observed=True)["ratio"].max().rename("zero_max_ratio")
    systems = sorted(df["system_id"].unique())
    sys_index = {s: i for i, s in enumerate(systems)}
    max_ratios = np.zeros((m, len(systems)), dtype=float)
    endpoint_sum: Dict[Tuple[str, str, str], float] = defaultdict(float)
    endpoint_exceed: Dict[Tuple[str, str, str], int] = defaultdict(int)
    endpoint_n: Dict[Tuple[str, str, str], int] = defaultdict(int)

    for rep in range(m):
        conc = impute_censored_values(df, fits, rng)
        loc_mean = aggregate_records(df, conc)
        ep = endpoint_ratios_from_location_means(loc_mean, benchmarks)
        sys_rat = ep.groupby("system_id", observed=True)["ratio"].max()
        for s, val in sys_rat.items():
            max_ratios[rep, sys_index[s]] = float(val)
        for row in ep.itertuples(index=False):
            key = (str(row.system_id), str(row.location_id), str(row.endpoint))
            val = float(row.ratio)
            endpoint_sum[key] += val
            endpoint_exceed[key] += int(val >= 1.0)
            endpoint_n[key] += 1
        if (rep + 1) % max(1, min(50, m // 10 if m >= 10 else 1)) == 0:
            print(f"  MI replicate {rep + 1}/{m}")

    mean_max = max_ratios.mean(axis=0)
    med_max = np.median(max_ratios, axis=0)
    p05 = np.percentile(max_ratios, 5, axis=0)
    p95 = np.percentile(max_ratios, 95, axis=0)
    p_exceed = (max_ratios >= 1.0).mean(axis=0)
    rows = []
    for s in systems:
        i = sys_index[s]
        z = float(zero_sys.get(s, 0.0))
        if flag_rule == "probability":
            mi_flag = p_exceed[i] >= prob_threshold
        else:
            mi_flag = mean_max[i] >= 1.0
        zero_flag = z >= 1.0
        if zero_flag and mi_flag:
            cls = "always_flagged"
        elif (not zero_flag) and mi_flag:
            cls = "mi_only_screening_flag"
        elif zero_flag and (not mi_flag):
            cls = "zero_only_flag_unexpected"
        else:
            cls = "never_flagged"
        rows.append({
            "system_id": s,
            "zero_max_ratio": z,
            "mi_mean_max_ratio": float(mean_max[i]),
            "mi_median_max_ratio": float(med_max[i]),
            "mi_p05_max_ratio": float(p05[i]),
            "mi_p95_max_ratio": float(p95[i]),
            "mi_exceedance_probability": float(p_exceed[i]),
            "delta_mean_minus_zero": float(mean_max[i] - z),
            "zero_flag": bool(zero_flag),
            "mi_flag": bool(mi_flag),
            "classification": cls,
        })
    systems_df = pd.DataFrame(rows).sort_values(["classification", "delta_mean_minus_zero"], ascending=[True, False])

    endpoint_rows = []
    zero_endpoint_keyed = zero_endpoint.set_index(["system_id", "location_id", "endpoint"])["ratio"].to_dict()
    for key, n in endpoint_n.items():
        mean_ratio = endpoint_sum[key] / n
        z = float(zero_endpoint_keyed.get(key, 0.0))
        endpoint_rows.append({
            "system_id": key[0],
            "location_id": key[1],
            "endpoint": key[2],
            "zero_ratio": z,
            "mi_mean_ratio": mean_ratio,
            "mi_exceedance_probability": endpoint_exceed[key] / n,
            "delta_mean_minus_zero": mean_ratio - z,
        })
    endpoints_df = pd.DataFrame(endpoint_rows).sort_values("delta_mean_minus_zero", ascending=False)

    if checkpoints is None:
        checkpoints = [50, 100, 200, 300, 500, 1000, m]
    conv_rows = []
    for c in sorted(set(int(x) for x in checkpoints if 1 <= int(x) <= m)):
        mm = max_ratios[:c, :].mean(axis=0)
        pe = (max_ratios[:c, :] >= 1.0).mean(axis=0)
        if flag_rule == "probability":
            mi_flags = pe >= prob_threshold
        else:
            mi_flags = mm >= 1.0
        zero_flags = np.array([float(zero_sys.get(s, 0.0)) >= 1.0 for s in systems])
        conv_rows.append({
            "replicates": c,
            "always_flagged": int(np.sum(zero_flags & mi_flags)),
            "mi_only_screening_flag": int(np.sum((~zero_flags) & mi_flags)),
            "never_flagged": int(np.sum((~zero_flags) & (~mi_flags))),
            "mean_mi_max_ratio": float(np.mean(mm)),
        })
    conv_df = pd.DataFrame(conv_rows)
    return {"systems": systems_df, "endpoints": endpoints_df, "convergence": conv_df, "zero_endpoints": zero_endpoint}


def fits_to_df(fits: Dict[str, FitResult]) -> pd.DataFrame:
    rows = []
    for f in fits.values():
        rows.append({
            "analyte": f.analyte,
            "tier": f.tier,
            "mode": f.mode,
            "n_total": f.n_total,
            "n_detect": f.n_detect,
            "n_censored": f.n_censored,
            "detection_rate_pct": 100 * f.n_detect / max(f.n_total, 1),
            "mu": f.mu,
            "sigma": f.sigma,
            "geometric_mean_ug_l": math.exp(f.mu),
            "arithmetic_mean_ug_l": math.exp(f.mu + 0.5 * f.sigma * f.sigma),
            "loglik": f.loglik,
            "success": f.success,
            "warnings": "; ".join(f.warnings),
            "eb_mu_bar": f.eb_mu_bar,
            "eb_s2": f.eb_s2,
            "eb_sigma_bar": f.eb_sigma_bar,
        })
    return pd.DataFrame(rows).sort_values(["tier", "analyte"])


def benchmark_to_df(benchmarks: Dict[str, Dict[str, float]]) -> pd.DataFrame:
    rows = []
    for k, v in benchmarks.get("individual_mcls_ug_l", {}).items():
        rows.append({"endpoint_type": "individual", "analyte": k, "benchmark_ug_l": v, "benchmark_ng_l": v * 1000})
    for k, v in benchmarks.get("hi_denominators_ug_l", {}).items():
        rows.append({"endpoint_type": "hazard_index_denominator", "analyte": k, "benchmark_ug_l": v, "benchmark_ng_l": v * 1000})
    return pd.DataFrame(rows)


def write_outputs(
    out_xlsx: Path,
    df: pd.DataFrame,
    qa_df: pd.DataFrame,
    fits: Dict[str, FitResult],
    benchmarks: Dict[str, Dict[str, float]],
    results: Dict[str, pd.DataFrame],
    config: Dict[str, Any],
) -> None:
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    fits_df = fits_to_df(fits)
    benchmarks_df = benchmark_to_df(benchmarks)
    systems_df = results["systems"]
    endpoints_df = results["endpoints"]
    conv_df = results["convergence"]
    summary = pd.DataFrame([
        {"metric": "tool_version", "value": TOOL_VERSION},
        {"metric": "run_utc", "value": datetime.now(timezone.utc).isoformat()},
        {"metric": "records_valid", "value": int(len(df))},
        {"metric": "systems", "value": int(df["system_id"].nunique())},
        {"metric": "locations", "value": int(df[["system_id", "location_id"]].drop_duplicates().shape[0])},
        {"metric": "analytes", "value": int(df["analyte"].nunique())},
        {"metric": "mi_replicates", "value": int(config["m"])},
        {"metric": "seed", "value": int(config["seed"])},
        {"metric": "always_flagged", "value": int((systems_df["classification"] == "always_flagged").sum())},
        {"metric": "mi_only_screening_flag", "value": int((systems_df["classification"] == "mi_only_screening_flag").sum())},
        {"metric": "never_flagged", "value": int((systems_df["classification"] == "never_flagged").sum())},
        {"metric": "max_delta_mean_minus_zero", "value": float(systems_df["delta_mean_minus_zero"].max()) if len(systems_df) else np.nan},
    ])
    config_df = pd.DataFrame([{"parameter": k, "value": json.dumps(v) if isinstance(v, (dict, list)) else v} for k, v in config.items()])
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        systems_df.to_excel(writer, sheet_name="Systems", index=False)
        endpoints_df.to_excel(writer, sheet_name="Endpoint_Details", index=False)
        fits_df.to_excel(writer, sheet_name="Fits", index=False)
        benchmarks_df.to_excel(writer, sheet_name="Benchmarks", index=False)
        qa_df.to_excel(writer, sheet_name="Import_QA", index=False)
        conv_df.to_excel(writer, sheet_name="Convergence", index=False)
        config_df.to_excel(writer, sheet_name="Config", index=False)
        for ws in writer.book.worksheets:
            ws.freeze_panes = "A2"
            for col in ws.columns:
                max_len = 0
                letter = col[0].column_letter
                for cell in col[:200]:
                    val = cell.value
                    if val is not None:
                        max_len = max(max_len, len(str(val)))
                ws.column_dimensions[letter].width = min(max(max_len + 2, 10), 45)
    stem = out_xlsx.with_suffix("")
    systems_df.to_csv(stem.with_name(stem.name + "_systems.csv"), index=False)
    endpoints_df.to_csv(stem.with_name(stem.name + "_endpoint_details.csv"), index=False)
    fits_df.to_csv(stem.with_name(stem.name + "_fits.csv"), index=False)
    with stem.with_name(stem.name + "_summary.json").open("w", encoding="utf-8") as f:
        json.dump({
            "summary": {r["metric"]: r["value"] for r in summary.to_dict(orient="records")},
            "config": config,
            "tier_counts": fits_df["tier"].value_counts().to_dict(),
            "benchmark_note": "Default benchmarks reflect the April 2024 PFAS NPDWR framework as configurable screening thresholds; verify current regulatory status before decision support.",
        }, f, indent=2)


def make_demo_data(path: Path, n_systems: int = 8, seed: int = 123) -> Path:
    rng = np.random.default_rng(seed)
    analytes = ["PFOA", "PFOS", "PFHXS", "PFNA", "HFPO-DA", "PFBS"]
    rows = []
    for s in range(n_systems):
        sysid = f"DEMO{str(s+1).zfill(4)}"
        n_locs = int(rng.integers(1, 4))
        sys_factor = rng.lognormal(0, 0.8)
        for loc in range(n_locs):
            for event in range(2):
                for a in analytes:
                    base = {"PFOA": 0.002, "PFOS": 0.0018, "PFHXS": 0.0015, "PFNA": 0.0012, "HFPO-DA": 0.001, "PFBS": 0.006}[a]
                    true = rng.lognormal(math.log(base * sys_factor), 0.9)
                    mrl = {"PFBS": 0.003}.get(a, 0.0015)
                    detect = true >= mrl
                    rows.append({
                        "PWSID": sysid,
                        "SamplingLocationID": f"{sysid}_EP{loc+1}",
                        "Analyte": a,
                        "Result": true if detect else f"<{mrl}",
                        "Detect": "Y" if detect else "N",
                        "MRL": mrl,
                    })
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="PFAS Censoring-Aware Screening Tool (PFAS-CAST), Python screening engine v1.1.1")
    p.add_argument("--input", type=Path, help="Input CSV. Omit when using --demo.")
    p.add_argument("--output", type=Path, default=Path("pfas_cast_screening_results_v1_1.xlsx"), help="Output .xlsx workbook.")
    p.add_argument("--unit", default="ug/L", help="Input concentration unit: ug/L or ng/L/ppt. Internal calculations use ug/L.")
    p.add_argument("--m", type=int, default=500, help="Number of MI replicates.")
    p.add_argument("--seed", type=int, default=42, help="Random seed.")
    p.add_argument("--benchmarks", type=Path, help="Optional JSON benchmark file with individual_mcls_ug_l and hi_denominators_ug_l.")
    p.add_argument("--flag-rule", choices=["mean", "probability"], default="mean", help="System flag rule for MI screening.")
    p.add_argument("--prob-threshold", type=float, default=0.5, help="Exceedance-probability threshold if --flag-rule=probability.")
    p.add_argument("--b2-boot", type=int, default=200, help="Conditional bootstrap draws for Tier B2 analytes.")
    p.add_argument("--b1-grid", type=int, default=40, help="Grid size per dimension for Tier B1 likelihood-weighted sampler.")
    p.add_argument("--demo", action="store_true", help="Generate and run a small synthetic demo dataset.")
    p.add_argument("--demo-csv", type=Path, default=Path("pfas_cast_demo.csv"), help="Where to write demo CSV when --demo is used.")
    p.add_argument("--system-col", dest="system_id", help="Explicit system ID column.")
    p.add_argument("--location-col", dest="location_id", help="Explicit location ID column.")
    p.add_argument("--analyte-col", dest="analyte", help="Explicit analyte column.")
    p.add_argument("--value-col", dest="value", help="Explicit result/concentration column.")
    p.add_argument("--detect-col", dest="detect", help="Explicit detect-flag column.")
    p.add_argument("--rl-col", dest="reporting_limit", help="Explicit reporting-limit/MRL column.")
    args = p.parse_args(argv)

    if args.demo:
        input_path = make_demo_data(args.demo_csv, seed=args.seed)
        print(f"Generated demo CSV: {input_path}")
    elif args.input:
        input_path = args.input
    else:
        p.error("Provide --input CSV or use --demo.")
        return 2
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    t0 = time.time()
    explicit_cols = {k: getattr(args, k) for k in COLUMN_ALIASES}
    df, qa_df, qa = parse_input_csv(input_path, args.unit, explicit_cols)
    benchmarks = load_benchmarks(args.benchmarks)
    print(f"Loaded {len(df):,} valid records from {input_path}")
    print(f"Systems={df['system_id'].nunique():,}; analytes={df['analyte'].nunique():,}; detections={int(df['detect'].sum()):,}")
    print("Fitting censored lognormal marginals...")
    fits = fit_all_analytes(df, seed=args.seed, b2_boot=args.b2_boot, b1_grid=args.b1_grid)
    tier_counts = pd.Series([f.tier for f in fits.values()]).value_counts().to_dict()
    print(f"Tier counts: {tier_counts}")
    print(f"Running record-level MI screening with M={args.m}...")
    results = run_screening(df, fits, benchmarks, args.m, args.seed, args.flag_rule, args.prob_threshold)
    config = {
        "tool_version": TOOL_VERSION,
        "input_file": str(input_path),
        "input_sha256": sha256_file(input_path),
        "unit_input": args.unit,
        "unit_internal": INTERNAL_UNIT,
        "m": args.m,
        "seed": args.seed,
        "flag_rule": args.flag_rule,
        "prob_threshold": args.prob_threshold,
        "benchmarks": benchmarks,
        "b2_boot": args.b2_boot,
        "b1_grid": args.b1_grid,
        "runtime_note": "Screening tool; not a regulatory-compliance determination.",
    }
    write_outputs(args.output, df, qa_df, fits, benchmarks, results, config)
    systems = results["systems"]
    print("\nClassification summary:")
    print(systems["classification"].value_counts().to_string())
    print(f"\nWrote: {args.output}")
    print(f"Runtime: {(time.time() - t0):.1f} seconds")
    print("\nRegulatory note: default benchmarks are configurable screening thresholds based on the April 2024 PFAS NPDWR framework; verify current status before decision support.")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
