"""Design-stage precision planning (NO patient data): Monte-Carlo SD of the interaction estimator for hypothetical
cohort sizes, drawing from the synthetic generator (median OS roughly 15 months wild-type, 9 months altered on
monotherapy; ~22% altered). Output is a planning aid, not a result.

    python -m geniebpc.precision outputs/diagnostics/design_precision.csv
"""
import sys
import numpy as np, pandas as pd
from . import synth, derive, estimators as E

SCENARIOS = {"A: 150 mono + 40 chemo-IO (pessimistic public-release guess)": (150, 40),
             "B: 250 mono + 100 chemo-IO": (250, 100),
             "C: 400 mono + 400 chemo-IO (larger consortium-scale guess)": (400, 400),
             "D: 1000 mono + 1000 chemo-IO": (1000, 1000)}


def main(out, reps=300):
    pool, _ = derive.build_cohort(synth.generate(n=120000, seed=99)[0])
    mono, chemo = pool[pool.treat == 0], pool[pool.treat == 1]
    rows = []
    for name, (n0, n1) in SCENARIOS.items():
        est = []
        for r in range(reps):
            d = pd.concat([mono.sample(n0, random_state=r), chemo.sample(n1, random_state=r)]).reset_index(drop=True)
            try:
                e = E.point_estimates(d); est.append((e["int_rmst12"], e["int_rmst24"], (d.biomarker[d.treat == 1]).sum()))
            except Exception:
                pass
        a = np.array(est, float); a = a[~np.isnan(a).any(axis=1)]
        for j, tau in enumerate((12, 24)):
            sd = a[:, j].std(ddof=1)
            rows.append(dict(scenario=name, horizon_months=tau, mc_sd_interaction_months=round(sd, 2),
                             ci_half_width_months=round(1.96 * sd, 2), detectable_interaction_80pct_power_months=round(2.8 * sd, 2),
                             mean_altered_in_chemo_arm=round(a[:, 2].mean(), 1), usable_reps=len(a)))
    pd.DataFrame(rows).to_csv(out, index=False); print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1])
