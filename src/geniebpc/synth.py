"""Synthetic GENIE-BPC-shaped data with a KNOWN treatment x biomarker effect. For pipeline validation only.
Nothing produced from these data is a study result."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from . import config as C

MONO = ["Pembrolizumab", "Nivolumab", "Atezolizumab"]
CHEMO = ["Carboplatin, Pembrolizumab, Pemetrexed Disodium", "Carboplatin, Paclitaxel, Pembrolizumab",
         "Atezolizumab, Bevacizumab, Carboplatin, Paclitaxel"]
OTHER = ["Carboplatin, Pemetrexed Disodium", "Erlotinib Hydrochloride", "Investigational Drug", "Ipilimumab, Nivolumab",
         "Carboplatin, Paclitaxel"]


def generate(n: int = 4000, seed: int = 1, b_treat: float = -0.25, b_bio: float = 0.55, b_int: float = -0.45,
             late_ngs_frac: float = 0.25, confounding: float = 1.0, max_fu_m: float = 48.0):
    """Returns (tables, truth). Hazard = h0*exp(b_bio*M + b_treat*A + b_int*A*M + covariates).
    truth holds the true RMST contrasts obtained from both potential outcomes of every ICI-treated patient."""
    r = np.random.default_rng(seed)
    rid = np.array([f"GENIE-SYN-{i:06d}" for i in range(n)])
    inst = r.choice(["MSK", "DFCI", "VICC", "UHN"], n, p=[.45, .35, .12, .08])
    age = np.clip(r.normal(66, 9, n), 35, 89)
    female = r.binomial(1, .52, n)
    adeno = r.binomial(1, .75, n)
    kras = r.binomial(1, np.where(adeno == 1, .35, .05))
    stk = r.binomial(1, np.where(kras == 1, .28, .10) * np.where(adeno == 1, 1, .3))
    keap = r.binomial(1, np.where(stk == 1, .45, .10))
    M = ((stk + keap) > 0).astype(int)
    pdl1 = np.where(r.random(n) < .2, np.nan, np.round(100 * r.beta(.5 - .15 * M, 1.2, n)))
    brain = r.binomial(1, .22, n); liver = r.binomial(1, .15, n)
    smoke = r.choice(["Never used", "Former user (quit >1 year)", "Current user"], n, p=[.15, .6, .25])
    de_novo = r.binomial(1, .7, n)
    ici = r.random(n) < .6
    high = np.nan_to_num(pdl1, nan=20) >= 50
    lp = confounding * (-1.6 * high + .9 * liver - .03 * (age - 66) + .3 * M) + .1
    A = r.binomial(1, 1 / (1 + np.exp(-lp)))
    xb = .5 * liver + .35 * brain + .02 * (age - 66) - .4 * high + .15 * (1 - adeno)
    u = r.exponential(1, n)

    def T(a):
        return u / (0.045 * np.exp(b_bio * M + b_treat * a + b_int * a * M + xb))
    T0, T1 = T(0), T(1)
    t = np.where(A == 1, T1, T0)
    cens = r.uniform(6, max_fu_m, n)
    os_m = np.minimum(t, cens); death = (t <= cens).astype(int)
    # NGS report timing relative to t0 (days); late reports create delayed entry: only survivors are sequenced
    rep_rel = np.where(r.random(n) < late_ngs_frac, r.uniform(15, 365, n), -r.uniform(1, 120, n)).round()
    keep = os_m * C.DAYS_PER_MONTH > np.maximum(rep_rel, 0)
    truth = {}
    for tau in C.TAUS_MONTHS:
        for m in (0, 1):
            i = ici & (M == m)
            truth[f"drmst{tau:g}_m{m}"] = float(np.mean(np.minimum(T1[i], tau) - np.minimum(T0[i], tau)))
        truth[f"int_rmst{tau:g}"] = truth[f"drmst{tau:g}_m1"] - truth[f"drmst{tau:g}_m0"]

    dmets_day = np.where(de_novo == 1, np.nan, r.integers(120, 900, n))
    adv = np.where(de_novo == 1, 0, dmets_day)
    t0 = adv + r.integers(14, 60, n)
    regimen = np.where(ici, np.where(A == 1, r.choice(CHEMO, n), r.choice(MONO, n)), r.choice(OTHER, n))
    assay = np.where(inst == "MSK", "MSK-IMPACT468", np.where(inst == "DFCI", "DFCI-ONCOPANEL-3",
                     np.where(inst == "VICC", "VICC-01-T7", "UHN-48-V1")))
    sample = np.char.add(rid, "-S1")
    k = keep
    pt = pd.DataFrame(dict(record_id=rid, institution=inst, birth_year=(2017 - age).round().astype(int),
                           naaccr_sex_code=np.where(female == 1, "Female", "Male"),
                           naaccr_race_code_primary=r.choice(["White", "Black", "Chinese", "Unknown"], n, p=[.75, .1, .08, .07])))[k]
    dx = pd.DataFrame(dict(record_id=rid, ca_seq=0, age_dx=age - t0 / 365.25, dob_ca_dx_days=((age - t0 / 365.25) * 365.25).round(),
                           stage_dx_iv=np.where(de_novo == 1, "Stage IV", "Stage I-III"),
                           dmets_stage_i_iii=np.where(de_novo == 1, np.nan, 1), dx_to_dmets_days=dmets_day,
                           ca_hist_adeno_squamous=np.where(adeno == 1, "Adenocarcinoma", "Squamous cell"),
                           ca_lung_cigarette=smoke, dmets_brain=brain, dx_to_dmets_brain_days=np.where(brain == 1, adv, np.nan),
                           dmets_liver=liver, dx_to_dmets_liver_days=np.where(liver == 1, adv, np.nan)))[k]
    reg1 = pd.DataFrame(dict(record_id=rid, ca_seq=0, regimen_number=1, regimen_drugs=regimen,
                             drugs_ct_yn=np.where(regimen == "Investigational Drug", "Yes", "No"), dx_reg_start_int=t0,
                             os_g_status=death, tt_os_g_days=(os_m * C.DAYS_PER_MONTH).round()))[k]
    # a later second-line regimen that must never be picked as first line
    l2 = reg1.sample(frac=.3, random_state=seed).assign(regimen_number=2, regimen_drugs="Docetaxel")
    l2["dx_reg_start_int"] += 200
    path = pd.DataFrame(dict(record_id=rid, dx_path_proc_days=t0 - r.integers(5, 40, n), pdl1_perc=pdl1))[k]
    # post-baseline PD-L1 result with a different value: using it would be leakage
    post = path.sample(frac=.2, random_state=seed).assign(pdl1_perc=99.0)
    post["dx_path_proc_days"] += 300
    cpt = pd.DataFrame(dict(record_id=rid, ca_seq=0, cpt_genie_sample_id=sample, cpt_seq_assay_id=assay,
                            dx_cpt_rep_days=t0 + rep_rel))[k]
    mut = []
    for g, flag, cls in [("STK11", stk, "Nonsense_Mutation"), ("KEAP1", keap, "Missense_Mutation"), ("KRAS", kras, "Missense_Mutation")]:
        mut.append(pd.DataFrame(dict(Hugo_Symbol=g, Variant_Classification=cls,
                                     HGVSp_Short="p.G12C" if g == "KRAS" else "p.X1", Tumor_Sample_Barcode=sample[(flag == 1) & k])))
    egfr = k & (r.random(n) < .05)
    mut.append(pd.DataFrame(dict(Hugo_Symbol="EGFR", Variant_Classification="Missense_Mutation", HGVSp_Short="p.L858R",
                                 Tumor_Sample_Barcode=sample[egfr])))
    tables = dict(pt_char=pt, ca_dx_index=dx, ca_drugs=pd.concat([reg1, l2]), prissmm_pathology=pd.concat([path, post]),
                  cpt=cpt, mutations_extended=pd.concat(mut))
    return {n_: t_.reset_index(drop=True) for n_, t_ in tables.items()}, truth


FILES = dict(pt_char="patient_level_dataset.csv", ca_dx_index="cancer_level_dataset_index.csv",
             ca_drugs="regimen_cancer_level_dataset.csv", prissmm_pathology="pathology_report_level_dataset.csv",
             cpt="cancer_panel_test_level_dataset.csv", mutations_extended="data_mutations_extended.txt")


def write(tables: dict, out: str | Path) -> None:
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        f = out / FILES[name]
        df.to_csv(f, sep="\t" if f.suffix == ".txt" else ",", index=False)
