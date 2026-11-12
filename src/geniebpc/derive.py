"""Source-to-derived-variable mapping: raw GENIE BPC tables -> one row per eligible patient.

All day counts are relative to the diagnosis date of the index NSCLC (the release's `dx_*` convention).
Time zero (t0) = start of the first systemic regimen for advanced disease (`dx_reg_start_int`).
Anything used as a baseline covariate or biomarker must be dated <= t0 (checked in tests/test_chronology.py).
"""
from __future__ import annotations
import re
import numpy as np
import pandas as pd
from . import config as C


# ---------------------------------------------------------------- treatment
def _norm(drug: str) -> str:
    return drug.split("(")[0].strip().lower()


def _is(drug: str, names: list[str], anywhere: bool = False) -> bool:
    return any((n in drug) if anywhere else drug.startswith(n) for n in names)


def classify_regimen(regimen_drugs: str, trial: str = "No") -> str:
    """ICI_mono | ICI_chemo | other:<reason>. Investigational/masked or trial regimens are never assigned an arm."""
    if not isinstance(regimen_drugs, str) or not regimen_drugs.strip():
        return "other:missing"
    drugs = [_norm(x) for x in re.sub(r"\([^)]*\)", "", regimen_drugs).split(",") if x.strip()]
    if str(trial).strip().lower() == "yes" or any("investigational" in d for d in drugs):
        return "other:trial_or_investigational"
    ici = [d for d in drugs if _is(d, C.PD1_PDL1)]
    plat = [d for d in drugs if _is(d, C.PLATINUM)]
    partner = [d for d in drugs if _is(d, C.PARTNER_CHEMO, anywhere=True) and d not in plat]
    adj = [d for d in drugs if _is(d, C.ALLOWED_ADJUNCT)]
    rest = [d for d in drugs if d not in ici + plat + partner + adj]
    if not ici:
        return "other:no_pd1_pdl1"
    if len(ici) > 1 or rest:
        return "other:ici_plus_nonprotocol_agent"
    if not plat and not partner and not adj:
        return C.ARM_MONO
    if plat:
        return C.ARM_CHEMO
    return "other:ici_nonplatinum_combo"


# ---------------------------------------------------------------- genomics
def _panel_covers_both(assay: pd.Series, covering: list[str] | None = None) -> pd.Series:
    pref = tuple(p.upper() for p in (covering or C.PANEL_PREFIX_BOTH_GENES))
    return assay.astype(str).str.upper().str.startswith(pref)


def sample_alterations(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """One row per sequenced sample with alteration flags (protein-altering and loss-of-function definitions)."""
    cpt = tables["cpt"][["record_id", "ca_seq", "cpt_genie_sample_id", "cpt_seq_assay_id", "dx_cpt_rep_days"]].copy()
    cpt = cpt.rename(columns={"cpt_genie_sample_id": "sample"})
    mut = tables["mutations_extended"]
    mut = mut[mut["Variant_Classification"].isin(C.PROTEIN_ALTERING)]
    cna = tables.get("cna")
    deep = {}
    if cna is not None and "Hugo_Symbol" in cna.columns:
        c = cna.set_index("Hugo_Symbol")
        for g in C.BIOMARKER_GENES:
            if g in c.index:
                row = pd.to_numeric(c.loc[[g]].iloc[0], errors="coerce")
                deep[g] = set(row.index[row <= -2])
    for g in C.BIOMARKER_GENES + ["KRAS", "TP53"]:
        mg = mut[mut["Hugo_Symbol"] == g]
        any_s = set(mg["Tumor_Sample_Barcode"])
        lof_s = set(mg.loc[mg["Variant_Classification"].isin(C.LOSS_OF_FUNCTION), "Tumor_Sample_Barcode"])
        cpt[f"{g}_alt"] = cpt["sample"].isin(any_s | deep.get(g, set()))
        cpt[f"{g}_lof"] = cpt["sample"].isin(lof_s | deep.get(g, set()))
    egfr = mut[mut["Hugo_Symbol"] == "EGFR"]
    hg = egfr["HGVSp_Short"].astype(str)
    sens = egfr[hg.str.contains(r"L858R|G719|L861Q|S768I", regex=True)
                | egfr["Variant_Classification"].isin(["In_Frame_Del", "In_Frame_Ins"])]
    drv = set(sens["Tumor_Sample_Barcode"])
    fus = tables.get("fusions")
    if fus is not None and {"Hugo_Symbol", "Tumor_Sample_Barcode"} <= set(fus.columns):
        drv |= set(fus.loc[fus["Hugo_Symbol"].isin(C.DRIVER_FUSION_GENES), "Tumor_Sample_Barcode"])
    cpt["driver"] = cpt["sample"].isin(drv)
    cpt["covers_both"] = _panel_covers_both(cpt["cpt_seq_assay_id"])
    return cpt


def _patient_biomarkers(samples: pd.DataFrame, fl: pd.DataFrame, pretreatment_only: bool) -> pd.DataFrame:
    s = samples.merge(fl[["record_id", "ca_seq", "t0"]], on=["record_id", "ca_seq"])
    s["rep_rel"] = s["dx_cpt_rep_days"] - s["t0"]
    out = s.groupby("record_id").agg(first_report_rel=("rep_rel", "min"), driver_any=("driver", "max"))
    ev = s[s["covers_both"]]
    out["first_evaluable_report_rel"] = ev.groupby("record_id")["rep_rel"].min()
    if pretreatment_only:
        ev = ev[ev["rep_rel"] <= 0]
    flags = [f"{g}_{k}" for g in C.BIOMARKER_GENES + ["KRAS", "TP53"] for k in ("alt", "lof")]
    out = out.join(ev.groupby("record_id")[flags].max()).join(ev.groupby("record_id").size().rename("n_evaluable"))
    return out.reset_index()


# ---------------------------------------------------------------- cohort
def first_line(tables: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, list[tuple[str, int]]]:
    """First systemic regimen started on/after the advanced-disease date of the first index NSCLC."""
    flow = []
    dx = tables["ca_dx_index"].sort_values(["record_id", "ca_seq"]).drop_duplicates("record_id")
    flow.append(("Patients with an index NSCLC diagnosis", dx["record_id"].nunique()))
    de_novo = dx["stage_dx_iv"].astype(str).eq("Stage IV")
    recur = (~de_novo) & (pd.to_numeric(dx["dmets_stage_i_iii"], errors="coerce") == 1)
    dx = dx.assign(de_novo_iv=de_novo, adv_day=np.where(de_novo, 0.0, pd.to_numeric(dx["dx_to_dmets_days"], errors="coerce")))
    dx = dx[de_novo | recur]
    dx = dx[dx["adv_day"].notna()]
    flow.append(("Advanced disease (stage IV at diagnosis or later distant metastasis)", len(dx)))
    reg = tables["ca_drugs"].merge(dx[["record_id", "ca_seq", "adv_day", "de_novo_iv"]], on=["record_id", "ca_seq"])
    reg["dx_reg_start_int"] = pd.to_numeric(reg["dx_reg_start_int"], errors="coerce")
    prior = reg[reg["dx_reg_start_int"] < reg["adv_day"]].groupby("record_id").size().rename("n_prior_regimens")
    reg = reg[reg["dx_reg_start_int"] >= reg["adv_day"]].sort_values(["record_id", "dx_reg_start_int", "regimen_number"])
    fl = reg.drop_duplicates("record_id").copy()
    flow.append(("Received >=1 systemic regimen for advanced disease", len(fl)))
    fl["arm_raw"] = [classify_regimen(d, t) for d, t in zip(fl["regimen_drugs"], fl["drugs_ct_yn"])]
    fl["t0"] = fl["dx_reg_start_int"]
    fl = fl.merge(prior.reset_index(), on="record_id", how="left").fillna({"n_prior_regimens": 0})
    return fl, flow


def _pdl1(tables, fl):
    p = tables["prissmm_pathology"][["record_id", "dx_path_proc_days", "pdl1_perc"]].copy()
    p["pdl1_perc"] = pd.to_numeric(p["pdl1_perc"], errors="coerce")
    p = p.dropna(subset=["pdl1_perc"]).merge(fl[["record_id", "t0"]], on="record_id")
    p = p[p["dx_path_proc_days"] <= p["t0"]].sort_values(["record_id", "dx_path_proc_days"])
    last = p.drop_duplicates("record_id", keep="last").set_index("record_id")
    cat = pd.cut(last["pdl1_perc"], [-0.1, 0.999, 49.999, 100], labels=["<1%", "1-49%", ">=50%"]).astype(str)
    return cat.rename("pdl1_cat"), last["dx_path_proc_days"].rename("pdl1_day")


def build_cohort(tables: dict[str, pd.DataFrame], pretreatment_only: bool = True,
                 lof_only: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (analytic dataset, patient-flow table).

    pretreatment_only=True (primary): biomarker status uses only NGS reports dated <= t0, so there is no
    immortal time between t0 and sequencing. False (sensitivity): any evaluable report; `entry_m` then carries
    the delayed-entry time (months from t0 to first evaluable report) for left-truncation-adjusted risk sets.
    """
    fl, flow = first_line(tables)
    fl = fl[fl["arm_raw"].isin([C.ARM_MONO, C.ARM_CHEMO])].copy()
    flow.append(("First-line regimen is anti-PD-(L)1 monotherapy or anti-PD-(L)1 + platinum chemotherapy "
                 "(non-trial, non-investigational)", len(fl)))
    bm = _patient_biomarkers(sample_alterations(tables), fl, pretreatment_only)
    d = fl.merge(bm, on="record_id", how="left")
    d = d[d["n_evaluable"].fillna(0) > 0]
    flow.append(("STK11 and KEAP1 evaluable on a covering NGS panel"
                 + (" reported on/before treatment start" if pretreatment_only else " (any time)"), len(d)))
    d = d[~d["driver_any"].fillna(False).astype(bool)]
    flow.append(("No sensitising EGFR mutation or ALK/ROS1 fusion detected", len(d)))
    d["os_m"] = pd.to_numeric(d["tt_os_g_days"], errors="coerce") / C.DAYS_PER_MONTH
    d["death"] = pd.to_numeric(d["os_g_status"], errors="coerce")
    d = d[d["os_m"].notna() & d["death"].isin([0, 1]) & (d["os_m"] > 0)]
    flow.append(("Valid overall-survival follow-up after time zero", len(d)))

    k = "lof" if lof_only else "alt"
    d["STK11"] = d[f"STK11_{k}"].astype(int)
    d["KEAP1"] = d[f"KEAP1_{k}"].astype(int)
    d["biomarker"] = ((d["STK11"] + d["KEAP1"]) > 0).astype(int)
    d["KRAS"] = d["KRAS_alt"].astype(int)
    d["TP53"] = d["TP53_alt"].astype(int)
    d["treat"] = (d["arm_raw"] == C.ARM_CHEMO).astype(int)
    d["entry_m"] = 0.0 if pretreatment_only else np.clip(d["first_evaluable_report_rel"], 0, None) / C.DAYS_PER_MONTH
    d = d[d["os_m"] > d["entry_m"]]
    if not pretreatment_only:
        flow.append(("Under observation after NGS report (delayed entry)", len(d)))

    dx = tables["ca_dx_index"].sort_values(["record_id", "ca_seq"]).drop_duplicates("record_id")
    pt = tables["pt_char"]
    d = d.merge(dx[["record_id", "age_dx", "ca_hist_adeno_squamous", "ca_lung_cigarette", "dmets_brain",
                    "dx_to_dmets_brain_days", "dmets_liver", "dx_to_dmets_liver_days", "dob_ca_dx_days"]], on="record_id")
    d = d.merge(pt[["record_id", "institution", "birth_year", "naaccr_sex_code", "naaccr_race_code_primary"]],
                on="record_id", suffixes=("_reg", ""))
    d["age"] = pd.to_numeric(d["age_dx"], errors="coerce") + d["t0"] / 365.25
    d["female"] = d["naaccr_sex_code"].astype(str).str.lower().eq("female").astype(int)
    race = d["naaccr_race_code_primary"].astype(str)
    d["race"] = np.select([race.eq("White"), race.eq("Black")], ["White", "Black"], "Other/unknown")
    h = d["ca_hist_adeno_squamous"].astype(str)
    d["histology"] = np.select([h.str.startswith("Adeno"), h.str.startswith("Squamous")], ["Adenocarcinoma", "Squamous"], "Other/NOS")
    s = d["ca_lung_cigarette"].astype(str)
    d["smoking"] = np.select([s.str.startswith("Never"), s.str.startswith("Current"), s.str.startswith("Former")],
                             ["Never", "Current", "Former"], "Unknown")
    for site in ("brain", "liver"):
        day = pd.to_numeric(d[f"dx_to_dmets_{site}_days"], errors="coerce")
        d[f"{site}_mets"] = ((pd.to_numeric(d[f"dmets_{site}"], errors="coerce") == 1) & (day <= d["t0"])).astype(int)
    d["de_novo"] = d["de_novo_iv"].astype(int)
    d["prior_systemic"] = (d["n_prior_regimens"] > 0).astype(int)
    # approximate calendar year of t0 (birth_year is year-only, so +-1 year error; see docs/limitations.md)
    d["year_t0"] = np.floor(pd.to_numeric(d["birth_year"], errors="coerce")
                            + (pd.to_numeric(d["dob_ca_dx_days"], errors="coerce") + d["t0"]) / 365.25)
    cat, day = _pdl1(tables, d)
    d = d.merge(cat.reset_index(), on="record_id", how="left").merge(day.reset_index(), on="record_id", how="left")
    d["pdl1_cat"] = d["pdl1_cat"].fillna("Unknown")
    keep = ["record_id", "institution", "treat", "biomarker", "STK11", "KEAP1", "KRAS", "TP53", "os_m", "death", "entry_m",
            "age", "female", "race", "histology", "smoking", "brain_mets", "liver_mets", "de_novo", "prior_systemic",
            "year_t0", "pdl1_cat", "pdl1_day", "t0", "first_report_rel", "first_evaluable_report_rel", "regimen_drugs"]
    flow_df = pd.DataFrame(flow, columns=["step", "n"])
    return d[keep].reset_index(drop=True), flow_df


COVARIATES_NUM = ["age", "year_t0"]
COVARIATES_BIN = ["female", "brain_mets", "liver_mets", "de_novo", "prior_systemic", "KRAS", "TP53"]
COVARIATES_CAT = ["institution", "race", "histology", "smoking", "pdl1_cat"]


def design_matrix(d: pd.DataFrame, include_biomarker: bool = True) -> pd.DataFrame:
    X = d[COVARIATES_NUM + COVARIATES_BIN].astype(float).copy()
    for c in COVARIATES_NUM:
        X[c] = X[c].fillna(X[c].median())
    X = pd.concat([X, pd.get_dummies(d[COVARIATES_CAT].astype(str), drop_first=True, dtype=float)], axis=1)
    if include_biomarker:
        X["biomarker"] = d["biomarker"].astype(float)
    X = X.loc[:, X.nunique() > 1]
    return X
