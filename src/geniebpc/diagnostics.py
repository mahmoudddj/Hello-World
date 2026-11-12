"""Balance, overlap and positivity diagnostics (aggregate outputs only)."""
from __future__ import annotations
import numpy as np
import pandas as pd
from .derive import design_matrix
from . import config as C


def _smd(x, a, w):
    def mv(g):
        ww = w[g]; m = np.average(x[g], weights=ww); return m, np.average((x[g] - m) ** 2, weights=ww)
    (m1, v1), (m0, v0) = mv(a == 1), mv(a == 0)
    s = np.sqrt((v1 + v0) / 2)
    return 0.0 if s == 0 else (m1 - m0) / s


def balance_table(d: pd.DataFrame, w: np.ndarray) -> pd.DataFrame:
    X = design_matrix(d, include_biomarker=False); a = d["treat"].values; rows = []
    for label, mask in [("overall", np.ones(len(d), bool)), ("biomarker_wt", (d["biomarker"] == 0).values),
                        ("biomarker_altered", (d["biomarker"] == 1).values)]:
        if mask.sum() < 4 or len(set(a[mask])) < 2:
            continue
        for c in X.columns:
            x = X[c].values[mask]
            rows.append(dict(stratum=label, covariate=c, smd_unweighted=_smd(x, a[mask], np.ones(mask.sum())),
                             smd_weighted=_smd(x, a[mask], w[mask])))
    return pd.DataFrame(rows)


def overlap_summary(d: pd.DataFrame, ps: np.ndarray, w: np.ndarray) -> pd.DataFrame:
    rows = []
    for m in (0, 1):
        for a in (0, 1):
            i = ((d["biomarker"] == m) & (d["treat"] == a)).values
            if i.sum() == 0:
                continue
            rows.append(dict(biomarker=m, treat=a, n=int(i.sum()), ps_min=ps[i].min(), ps_p05=np.quantile(ps[i], .05),
                             ps_median=np.median(ps[i]), ps_p95=np.quantile(ps[i], .95), ps_max=ps[i].max(),
                             max_weight=w[i].max(), ess=w[i].sum() ** 2 / (w[i] ** 2).sum()))
    return pd.DataFrame(rows)


def cell_counts(d: pd.DataFrame) -> pd.DataFrame:
    """Treatment x biomarker cells: n, deaths, person-months, deaths by 12/24 months; drives the G1/G2 gate."""
    g = d.assign(arm=np.where(d["treat"] == 1, C.ARM_CHEMO, C.ARM_MONO),
                 status=np.where(d["biomarker"] == 1, "STK11/KEAP1 altered", "STK11/KEAP1 wild-type"))
    t = g.groupby(["status", "arm"]).agg(n=("death", "size"), deaths=("death", "sum"),
                                         person_months=("os_m", "sum"), median_followup_m=("os_m", "median")).reset_index()
    for tau in C.TAUS_MONTHS:
        t[f"deaths_by_{tau:g}m"] = g.assign(e=(g["death"] == 1) & (g["os_m"] <= tau)).groupby(["status", "arm"])["e"].sum().values
        t[f"at_risk_beyond_{tau:g}m"] = g.assign(e=g["os_m"] > tau).groupby(["status", "arm"])["e"].sum().values
    return t


def gate_positivity(cells: pd.DataFrame) -> tuple[bool, str]:
    if len(cells) < 4:
        return False, "at least one treatment x biomarker cell is empty"
    bad = cells[(cells["n"] < C.MIN_N_PER_CELL) | (cells["deaths"] < C.MIN_EVENTS_PER_CELL)]
    if len(bad):
        return False, (f"{len(bad)} of 4 cells below pre-specified minimum "
                       f"(n>={C.MIN_N_PER_CELL}, deaths>={C.MIN_EVENTS_PER_CELL})")
    return True, "all four treatment x biomarker cells meet the pre-specified minimum"
