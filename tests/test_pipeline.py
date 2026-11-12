"""Automated checks on synthetic fixtures (no patient data). Run: pytest -q"""
import numpy as np, pandas as pd, pytest
from lifelines import KaplanMeierFitter
from lifelines.utils import restricted_mean_survival_time
from geniebpc import config as C, derive, estimators as E, diagnostics as D, synth, io


@pytest.fixture(scope="module")
def syn():
    return synth.generate(n=6000, seed=7)


# ---------------- treatment classification
@pytest.mark.parametrize("drugs,trial,expected", [
    ("Pembrolizumab", "No", C.ARM_MONO),
    ("Nivolumab(BMS936558,MDX1106,NIVO,ONO4538,Opdivo)", "No", C.ARM_MONO),
    ("Carboplatin, Pembrolizumab, Pemetrexed Disodium", "No", C.ARM_CHEMO),
    ("Carboplatin, Paclitaxel Protein-Bound, Pembrolizumab", "No", C.ARM_CHEMO),
    ("Atezolizumab, Bevacizumab, Carboplatin, Paclitaxel", "No", C.ARM_CHEMO),
    ("Ipilimumab, Nivolumab", "No", "other:ici_plus_nonprotocol_agent"),
    ("Carboplatin, Pemetrexed Disodium", "No", "other:no_pd1_pdl1"),
    ("Pembrolizumab, Pemetrexed Disodium", "No", "other:ici_nonplatinum_combo"),
    ("Investigational Drug, Pembrolizumab", "No", "other:trial_or_investigational"),
    ("Pembrolizumab", "Yes", "other:trial_or_investigational"),
    ("Erlotinib Hydrochloride, Nivolumab", "No", "other:ici_plus_nonprotocol_agent"),
    (None, "No", "other:missing"),
])
def test_classify_regimen(drugs, trial, expected):
    assert derive.classify_regimen(drugs, trial) == expected


# ---------------- estimators against an independent implementation
def test_km_and_rmst_match_lifelines():
    r = np.random.default_rng(0); t = r.exponential(10, 400).round(1) + .1; e = r.binomial(1, .7, 400)
    k = KaplanMeierFitter().fit(t, e)
    for tau in (6, 12, 24):
        assert E.rmst(t, e, tau) == pytest.approx(restricted_mean_survival_time(k, t=tau), abs=1e-6)
        assert E.surv_at(t, e, tau) == pytest.approx(float(k.predict(tau)), abs=1e-9)


def test_weighted_km_equals_replication():
    r = np.random.default_rng(1); t = r.exponential(10, 200).round(1) + .1; e = r.binomial(1, .7, 200); w = r.integers(1, 4, 200)
    assert E.rmst(t, e, 12, w) == pytest.approx(E.rmst(np.repeat(t, w), np.repeat(e, w), 12), abs=1e-9)


def test_delayed_entry_km_matches_lifelines():
    r = np.random.default_rng(2); t = r.exponential(10, 500) + .5; e = r.binomial(1, .8, 500); ent = np.minimum(r.uniform(0, 3, 500), t - .1)
    ent[:150] = 0
    k = KaplanMeierFitter().fit(t, e, entry=ent)
    ts, s = E.km(t, e, entry=ent)
    assert s[-1] == pytest.approx(float(k.survival_function_.iloc[-1, 0]), abs=1e-9)


def test_rmst_no_events_is_tau():
    assert E.rmst(np.array([5., 9.]), np.array([0, 0]), 4.0) == 4.0


# ---------------- chronology / leakage
def test_first_line_is_earliest_regimen_and_arms_are_protocol(syn):
    d, _ = derive.build_cohort(syn[0])
    assert d["regimen_drugs"].ne("Docetaxel").all()
    assert set(d["treat"]) == {0, 1} and d["record_id"].is_unique


def test_no_post_baseline_information_in_primary_cohort(syn):
    d, _ = derive.build_cohort(syn[0], pretreatment_only=True)
    assert (d["first_evaluable_report_rel"] <= 0).all()          # biomarker known at time zero
    assert (d["entry_m"] == 0).all()
    assert (d["pdl1_day"].dropna() <= d.loc[d["pdl1_day"].notna(), "t0"]).all()
    assert not (d["pdl1_cat"].eq(">=50%") & d["pdl1_day"].isna()).any()


def test_post_baseline_ngs_is_ignored_for_primary_biomarker():
    tabs, _ = synth.generate(n=800, seed=3, late_ngs_frac=1.0)       # every report after t0
    d, _ = derive.build_cohort(tabs, pretreatment_only=True)
    assert len(d) == 0
    d2, _ = derive.build_cohort(tabs, pretreatment_only=False)
    assert len(d2) > 0 and (d2["entry_m"] > 0).all() and (d2["os_m"] > d2["entry_m"]).all()


def test_exclusions(syn):
    tabs = syn[0]; d, flow = derive.build_cohort(tabs)
    egfr = set(tabs["mutations_extended"].query("Hugo_Symbol=='EGFR'")["Tumor_Sample_Barcode"].str[:-3])
    assert not (set(d["record_id"]) & egfr)
    assert d["institution"].ne("UHN").all()                          # hotspot panel: KEAP1 not evaluable
    assert flow["n"].is_monotonic_decreasing


def test_lof_definition_is_stricter(syn):
    a, _ = derive.build_cohort(syn[0]); b, _ = derive.build_cohort(syn[0], lof_only=True)
    assert b["KEAP1"].sum() == 0 and b["STK11"].sum() == a["STK11"].sum()   # synthetic KEAP1 are missense


# ---------------- causal estimators recover the planted truth
def test_iptw_and_aipw_recover_truth_and_improve_balance():
    tabs, truth = synth.generate(n=40000, seed=21); d, _ = derive.build_cohort(tabs)
    est = E.point_estimates(d)
    for k in truth:                                   # large sample: IPTW must be close to the planted truth
        assert est[k] == pytest.approx(truth[k], abs=0.45), k
    bal = D.balance_table(d, E.iptw(d)).query("stratum=='overall'")
    assert bal["smd_weighted"].abs().max() < 0.1 < bal["smd_unweighted"].abs().max()
    dr = E.aipw_estimates(d.sample(5000, random_state=0).reset_index(drop=True))   # smaller n: looser tolerance
    for k in truth:
        assert dr[k] == pytest.approx(truth[k], abs=1.3), k


def test_null_interaction_not_manufactured():
    tabs, truth = synth.generate(n=8000, seed=11, b_int=0.0, b_bio=0.0)
    d, _ = derive.build_cohort(tabs); est = E.point_estimates(d)
    assert abs(truth["int_rmst24"]) < 0.05 and abs(est["int_rmst24"]) < 1.5


def test_prognostic_only_biomarker_gives_no_hr_interaction():
    """Strong prognostic effect, no effect modification on the hazard scale -> product-term HR ~ 1."""
    tabs, _ = synth.generate(n=8000, seed=5, b_int=0.0, b_bio=0.8)
    d, _ = derive.build_cohort(tabs); c = E.cox_interaction(d, E.iptw(d)).set_index("term")
    assert c.loc["biomarker", "HR"] > 1.6 and 0.75 < c.loc["treat_x_biomarker", "HR"] < 1.33


def test_delayed_entry_removes_immortal_time_bias():
    tabs, _ = synth.generate(n=8000, seed=9, late_ngs_frac=0.9, b_int=0, b_bio=0, b_treat=0, confounding=0)
    d, _ = derive.build_cohort(tabs, pretreatment_only=False)
    naive = E.rmst(d["os_m"].values, d["death"].values, 12)
    adj = E.rmst(d["os_m"].values, d["death"].values, 12, entry=d["entry_m"].values)
    ref, _ = derive.build_cohort(synth.generate(n=8000, seed=9, late_ngs_frac=0.0, b_int=0, b_bio=0, b_treat=0, confounding=0)[0])
    true = E.rmst(ref["os_m"].values, ref["death"].values, 12)
    assert naive > true + 0.4 and abs(adj - true) < abs(naive - true)


# ---------------- gates, IO, reproducibility
def test_gate_fails_on_sparse_cells():
    tabs, _ = synth.generate(n=250, seed=2); d, _ = derive.build_cohort(tabs)
    ok, why = D.gate_positivity(D.cell_counts(d)); assert not ok and "cell" in why


def test_schema_report_flags_missing_column(tmp_path):
    tabs, _ = synth.generate(n=200, seed=1); tabs["ca_drugs"] = tabs["ca_drugs"].drop(columns="tt_os_g_days")
    synth.write(tabs, tmp_path); _, rep = io.load_tables(str(tmp_path))
    assert not io.schema_ok(rep) and "tt_os_g_days" in set(rep.loc[~rep.present, "column"])


def test_bootstrap_is_deterministic(syn):
    d, _ = derive.build_cohort(syn[0]); d = d.head(500)
    a = E.bootstrap(d, E.point_estimates, B=5); b = E.bootstrap(d, E.point_estimates, B=5)
    pd.testing.assert_frame_equal(a, b)


def test_evalue():
    assert E.e_value(1.0) == pytest.approx(1.0) and E.e_value(0.5) == pytest.approx(E.e_value(2.0)) and E.e_value(2.0) > 1.5
