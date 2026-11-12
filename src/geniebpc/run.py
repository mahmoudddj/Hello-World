"""Entry point: G1 feasibility -> gate -> (only if the gate passes) primary and sensitivity analyses.

    python -m geniebpc.run --data-dir "$GENIE_BPC_DATA_DIR" --out outputs [--boot 1000]
    python -m geniebpc.run --synthetic --out outputs_synthetic     # pipeline validation, not study results
"""
from __future__ import annotations
import argparse, json, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from . import config as C, io, derive, estimators as E, diagnostics as D, synth


def _md(df: pd.DataFrame) -> str:
    return df.to_markdown(index=False, floatfmt=".3g")


def feasibility(tables, schema, out: Path, label: str) -> tuple[pd.DataFrame | None, bool, str]:
    (out / "tables").mkdir(parents=True, exist_ok=True); (out / "diagnostics").mkdir(parents=True, exist_ok=True)
    schema.to_csv(out / "diagnostics" / "schema_check.csv", index=False)
    L = [f"# G1 feasibility report ({label})", "", "## Schema check", "",
         f"Required tables/columns present: **{io.schema_ok(schema)}**", "", _md(schema[~schema.present]) if (~schema.present).any() else "No missing items.", ""]
    if not io.schema_ok(schema):
        (out / "feasibility_auto.md").write_text("\n".join(L)); return None, False, "required tables or columns missing"
    fl, _ = derive.first_line(tables)
    arms = fl["arm_raw"].value_counts().rename_axis("first_line_class").reset_index(name="n")
    d, flow = derive.build_cohort(tables, pretreatment_only=True)
    d_any, flow_any = derive.build_cohort(tables, pretreatment_only=False)
    cells = D.cell_counts(d); ok, why = D.gate_positivity(cells)
    flow.to_csv(out / "tables" / "cohort_flow_primary.csv", index=False)
    flow_any.to_csv(out / "tables" / "cohort_flow_delayed_entry.csv", index=False)
    cells.to_csv(out / "tables" / "cells_treatment_by_biomarker.csv", index=False)
    miss = pd.DataFrame({"variable": derive.COVARIATES_NUM + derive.COVARIATES_BIN + derive.COVARIATES_CAT})
    miss["missing_or_unknown_pct"] = [100 * float((d[v].isna() | d[v].astype(str).isin(["Unknown", "nan"])).mean()) if len(d) else np.nan
                                      for v in miss["variable"]]
    chron = pd.DataFrame(dict(
        quantity=["NGS report on/before t0 (any-time cohort), %", "median days from t0 to first evaluable NGS report",
                  "PD-L1 result dated after t0 used as baseline (must be 0)", "year of t0: min", "year of t0: max"],
        value=[100 * float((d_any["first_evaluable_report_rel"] <= 0).mean()) if len(d_any) else np.nan,
               float(d_any["first_evaluable_report_rel"].median()) if len(d_any) else np.nan,
               int((d["pdl1_day"] > d["t0"]).sum()), d["year_t0"].min(), d["year_t0"].max()]))
    xt = pd.crosstab(d["STK11"].map({0: "STK11 wt", 1: "STK11 alt"}), d["KEAP1"].map({0: "KEAP1 wt", 1: "KEAP1 alt"})).reset_index()
    L += ["## First-line regimen classes (all advanced patients)", "", _md(arms), "", "## Cohort flow (primary: NGS report on/before time zero)", "",
          _md(flow), "", "## Cohort flow (sensitivity: any-time NGS with delayed entry)", "", _md(flow_any), "",
          "## Treatment x biomarker cells (primary cohort)", "", _md(cells), "", "## STK11 x KEAP1", "", _md(xt), "",
          "## Covariate missingness (primary cohort)", "", _md(miss), "", "## Chronology", "", _md(chron), "",
          f"## Positivity / precision gate\n\n**{'PASS' if ok else 'FAIL'}** - {why}.", ""]
    (out / "feasibility_auto.md").write_text("\n".join(L))
    return d, ok, why


def analyse(d: pd.DataFrame, tables, out: Path, B: int) -> None:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    (out / "figures").mkdir(parents=True, exist_ok=True)
    ps = E.propensity(d); w = E.iptw(d, ps)
    D.balance_table(d, w).to_csv(out / "diagnostics" / "balance_smd.csv", index=False)
    D.overlap_summary(d, ps, w).to_csv(out / "diagnostics" / "overlap_weights.csv", index=False)
    res = {"primary_iptw": E.bootstrap(d, E.point_estimates, B=B), "unadjusted": E.bootstrap(d, E.point_estimates, B=B, adjust=False),
           "aipw_pseudo": E.bootstrap(d, E.aipw_estimates, B=min(B, 200)),
           "sens_stratum_specific_ps": E.bootstrap(d, E.point_estimates, B=B, stratified_ps=True)}
    for name, kw in [("sens_delayed_entry_anytime_ngs", dict(pretreatment_only=False)), ("sens_loss_of_function_only", dict(lof_only=True))]:
        ds, _ = derive.build_cohort(tables, **kw)
        if D.gate_positivity(D.cell_counts(ds))[0]:
            res[name] = E.bootstrap(ds, E.point_estimates, B=B)
    for gene in ("STK11", "KEAP1"):
        dg = d.assign(biomarker=d[gene])
        if D.gate_positivity(D.cell_counts(dg))[0]:
            res[f"secondary_{gene}_alone"] = E.bootstrap(dg, E.point_estimates, B=B)
    pd.concat([r.assign(analysis=k) for k, r in res.items()]).to_csv(out / "tables" / "effect_estimates.csv", index=False)
    E.cox_interaction(d, w).to_csv(out / "tables" / "cox_interaction_secondary.csv", index=False)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for j, m in enumerate((0, 1)):
        for a, col in ((0, "#1f77b4"), (1, "#d62728")):
            i = ((d["biomarker"] == m) & (d["treat"] == a)).values
            t, s = E.km(d["os_m"].values[i], d["death"].values[i], w[i])
            ax[j].step(np.r_[0, t], np.r_[1, s], where="post", color=col, label=[C.ARM_MONO, C.ARM_CHEMO][a] + f" (n={i.sum()})")
        ax[j].set_title(["STK11/KEAP1 wild-type", "STK11 and/or KEAP1 altered"][m]); ax[j].set_xlabel("Months from first-line start")
        ax[j].set_xlim(0, 36); ax[j].legend(frameon=False)
    ax[0].set_ylabel("IPTW-adjusted overall survival"); fig.tight_layout(); fig.savefig(out / "figures" / "km_iptw_by_biomarker.png", dpi=200)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    for a, col in ((0, "#1f77b4"), (1, "#d62728")):
        ax.hist(ps[d["treat"].values == a], bins=25, alpha=.5, color=col, label=[C.ARM_MONO, C.ARM_CHEMO][a])
    ax.set_xlabel("Propensity of ICI + chemotherapy"); ax.legend(frameon=False); fig.tight_layout()
    fig.savefig(out / "figures" / "propensity_overlap.png", dpi=200)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir"); ap.add_argument("--out", default="outputs"); ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--synthetic", action="store_true"); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(argv); out = Path(a.out)
    truth = None
    if a.synthetic:
        tabs, truth = synth.generate(seed=a.seed)
        with tempfile.TemporaryDirectory() as tmp:
            synth.write(tabs, tmp); tables, schema = io.load_tables(tmp)
        label = "SYNTHETIC DATA - pipeline validation only, not study results"
    else:
        tables, schema = io.load_tables(a.data_dir); label = "authorised GENIE BPC data"
    d, ok, why = feasibility(tables, schema, out, label)
    status = dict(label=label, gate_pass=ok, reason=why, n_primary=None if d is None else int(len(d)), synthetic_truth=truth)
    if ok:
        analyse(d, tables, out, a.boot)
    else:
        print(f"GATE FAIL: {why}. Primary causal analysis NOT run (see {out}/feasibility_auto.md).")
    (out / "run_status.json").write_text(json.dumps(status, indent=2))
    print(json.dumps(status, indent=2)); return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
