"""Estimators: weighted Kaplan-Meier with delayed entry, RMST, IPTW, pseudo-observation AIPW, bootstrap.

Estimand (protocol/SAP.md): Delta(m, tau) = RMST_tau(ICI+chemo | M=m) - RMST_tau(ICI mono | M=m) in the
covariate-standardised population of biomarker stratum m; interaction = Delta(1, tau) - Delta(0, tau).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.preprocessing import StandardScaler
from .derive import design_matrix


def km(time, event, weights=None, entry=None):
    """Weighted Kaplan-Meier. Risk set at t: entry < t <= time. Returns (event_times, survival)."""
    time = np.asarray(time, float); event = np.asarray(event, float)
    w = np.ones_like(time) if weights is None else np.asarray(weights, float)
    o = np.argsort(time); ts, cw = time[o], np.cumsum(w[o])
    tot = cw[-1] if len(cw) else 0.0
    ut = np.unique(time[event == 1])
    if len(ut) == 0:
        return ut, np.ones(0)
    at_risk = tot - np.concatenate([[0.0], cw])[np.searchsorted(ts, ut, side="left")]
    if entry is not None:
        e = np.asarray(entry, float); oe = np.argsort(e); es, ce = e[oe], np.cumsum(w[oe])
        at_risk = at_risk - (ce[-1] - np.concatenate([[0.0], ce])[np.searchsorted(es, ut, side="left")])
    ev = event == 1
    oe2 = np.argsort(time[ev]); te, cd = time[ev][oe2], np.cumsum(w[ev][oe2])
    d = np.diff(np.concatenate([[0.0], cd[np.searchsorted(te, ut, side="right") - 1]]))
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(at_risk > 0, 1 - d / at_risk, 1.0)
    return ut, np.cumprod(np.clip(f, 0, 1))


def rmst(time, event, tau, weights=None, entry=None) -> float:
    t, s = km(time, event, weights, entry)
    keep = t < tau
    grid = np.concatenate([[0.0], t[keep], [tau]])
    return float(np.sum(np.diff(grid) * np.concatenate([[1.0], s[keep]])))


def surv_at(time, event, tau, weights=None, entry=None) -> float:
    t, s = km(time, event, weights, entry)
    return float(s[t <= tau][-1]) if np.any(t <= tau) else 1.0


def propensity(d: pd.DataFrame, C_reg: float = 10.0, stratified: bool = False, min_arm: int = 40) -> np.ndarray:
    """P(ICI+chemo | covariates, biomarker). stratified=True (pre-specified sensitivity) fits a separate model in
    each biomarker stratum when both arms have >= min_arm patients there, otherwise falls back to the pooled model."""
    if stratified:
        ps = propensity(d, C_reg)
        for m in (0, 1):
            i = (d["biomarker"] == m).values
            if min(d.loc[i, "treat"].sum(), (1 - d.loc[i, "treat"]).sum()) >= min_arm:
                ps[i] = propensity(d[i].assign(biomarker=d["biomarker"].values[i]), C_reg)
        return ps
    X = design_matrix(d)
    Z = StandardScaler().fit_transform(X)
    m = LogisticRegression(C=C_reg, max_iter=2000).fit(Z, d["treat"].values)
    return m.predict_proba(Z)[:, 1]


def iptw(d: pd.DataFrame, ps: np.ndarray | None = None, trunc=(0.01, 0.99)) -> np.ndarray:
    """Stabilised ATE weights, stabilised within biomarker stratum, truncated at given weight percentiles."""
    ps = propensity(d) if ps is None else ps
    a, m = d["treat"].values, d["biomarker"].values
    pa = pd.Series(a).groupby(m).transform("mean").values
    w = np.where(a == 1, pa / ps, (1 - pa) / (1 - ps))
    lo, hi = np.quantile(w, trunc)
    return np.clip(w, lo, hi)


def _contrasts(d, w, taus, use_entry):
    out = {}
    ent = d["entry_m"].values if use_entry else None
    for tau in taus:
        for m in (0, 1):
            for a in (0, 1):
                i = ((d["biomarker"] == m) & (d["treat"] == a)).values
                e = None if ent is None else ent[i]
                wi = None if w is None else w[i]
                if i.sum() == 0:
                    out[f"rmst{tau:g}_m{m}_a{a}"] = np.nan; out[f"surv{tau:g}_m{m}_a{a}"] = np.nan; continue
                out[f"rmst{tau:g}_m{m}_a{a}"] = rmst(d["os_m"].values[i], d["death"].values[i], tau, wi, e)
                out[f"surv{tau:g}_m{m}_a{a}"] = surv_at(d["os_m"].values[i], d["death"].values[i], tau, wi, e)
            for k in ("rmst", "surv"):
                out[f"d{k}{tau:g}_m{m}"] = out[f"{k}{tau:g}_m{m}_a1"] - out[f"{k}{tau:g}_m{m}_a0"]
        for k in ("rmst", "surv"):
            out[f"int_{k}{tau:g}"] = out[f"d{k}{tau:g}_m1"] - out[f"d{k}{tau:g}_m0"]
    return out


def point_estimates(d: pd.DataFrame, taus=(12.0, 24.0), adjust: bool = True, stratified_ps: bool = False) -> dict:
    use_entry = bool((d["entry_m"] > 0).any())
    w = iptw(d, propensity(d, stratified=stratified_ps)) if adjust else None
    return _contrasts(d, w, taus, use_entry)


# ---------------------------------------------------------------- doubly robust (pseudo-observations)
def pseudo_obs(time, event, tau) -> np.ndarray:
    """Leave-one-out jackknife pseudo-observations of RMST(tau) (no delayed entry)."""
    n = len(time); full = rmst(time, event, tau)
    idx = np.arange(n); loo = np.empty(n)
    for i in range(n):
        k = idx != i
        loo[i] = rmst(time[k], event[k], tau)
    return n * full - (n - 1) * loo


def aipw_estimates(d: pd.DataFrame, taus=(12.0, 24.0)) -> dict:
    """AIPW on RMST pseudo-observations (computed within arm x biomarker so that censoring may depend on both).
    Consistent if either the propensity model or the outcome regression is right. Not valid with delayed entry."""
    if (d["entry_m"] > 0).any():
        raise ValueError("pseudo-observation AIPW is only defined for the no-delayed-entry (primary) cohort")
    X = design_matrix(d).values
    a, m = d["treat"].values, d["biomarker"].values
    ps = np.clip(propensity(d), 0.02, 0.98)
    out = {}
    for tau in taus:
        theta = np.full(len(d), np.nan)
        for mm in (0, 1):
            for aa in (0, 1):
                i = (m == mm) & (a == aa)
                if i.sum() > 1:
                    theta[i] = pseudo_obs(d["os_m"].values[i], d["death"].values[i], tau)
        mu = {}
        for aa in (0, 1):
            i = (a == aa) & ~np.isnan(theta)
            mu[aa] = LinearRegression().fit(X[i], theta[i]).predict(X)
        phi1 = mu[1] + (a == 1) / ps * (np.nan_to_num(theta) - mu[1])
        phi0 = mu[0] + (a == 0) / (1 - ps) * (np.nan_to_num(theta) - mu[0])
        for mm in (0, 1):
            out[f"drmst{tau:g}_m{mm}"] = float(np.mean((phi1 - phi0)[m == mm]))
        out[f"int_rmst{tau:g}"] = out[f"drmst{tau:g}_m1"] - out[f"drmst{tau:g}_m0"]
    return out


# ---------------------------------------------------------------- bootstrap
def bootstrap(d: pd.DataFrame, fn, B: int = 500, seed: int = 20261112, **kw) -> pd.DataFrame:
    """Nonparametric patient-level bootstrap; every replicate re-fits the nuisance models. Returns a table with
    estimate, bootstrap SE, percentile 95% CI and a Wald p-value (reported for interaction terms only)."""
    from scipy.stats import norm
    rng = np.random.default_rng(seed)
    est = fn(d, **kw)
    reps = []
    for _ in range(B):
        b = d.iloc[rng.integers(0, len(d), len(d))].reset_index(drop=True)
        try:
            reps.append(fn(b, **kw))
        except Exception:
            continue
    R = pd.DataFrame(reps)
    rows = []
    for k, v in est.items():
        r = R[k].dropna() if k in R else pd.Series(dtype=float)
        se = r.std(ddof=1) if len(r) > 1 else np.nan
        lo, hi = (np.quantile(r, [0.025, 0.975]) if len(r) > 1 else (np.nan, np.nan))
        p = 2 * norm.sf(abs(v / se)) if (k.startswith(("int_", "d")) and se and se > 0) else np.nan
        rows.append(dict(parameter=k, estimate=v, se=se, ci_low=lo, ci_high=hi, p_wald=p, n_boot=len(r)))
    return pd.DataFrame(rows)


def cox_interaction(d: pd.DataFrame, weights: np.ndarray | None = None) -> pd.DataFrame:
    """Secondary: weighted Cox model with treat x biomarker product term and robust variance."""
    from lifelines import CoxPHFitter
    x = pd.DataFrame(dict(os_m=d["os_m"], death=d["death"], treat=d["treat"], biomarker=d["biomarker"],
                          treat_x_biomarker=d["treat"] * d["biomarker"], entry_m=d["entry_m"]))
    x["w"] = 1.0 if weights is None else weights
    delayed = bool((x["entry_m"] > 0).any())
    if not delayed:
        x = x.drop(columns="entry_m")
    f = CoxPHFitter().fit(x, "os_m", "death", weights_col="w", robust=True, entry_col="entry_m" if delayed else None)
    s = f.summary[["coef", "se(coef)", "coef lower 95%", "coef upper 95%", "p"]].copy()
    s["HR"] = np.exp(s["coef"]); s["HR_low"] = np.exp(s["coef lower 95%"]); s["HR_high"] = np.exp(s["coef upper 95%"])
    return s.reset_index().rename(columns={"covariate": "term"})


def e_value(hr: float) -> float:
    """E-value for a hazard ratio with a common outcome (VanderWeele & Ding 2017 approximation via RR ~ HR-based)."""
    rr = (1 - 0.5 ** np.sqrt(hr)) / (1 - 0.5 ** np.sqrt(1 / hr)) if hr != 1 else 1.0
    rr = 1 / rr if rr < 1 else rr
    return float(rr + np.sqrt(rr * (rr - 1)))
