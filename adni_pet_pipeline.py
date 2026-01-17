import pandas as pd
import numpy as np

# -----------------------
# Utilities
# -----------------------
def load_csv(path, upper_cols=True, parse_dates=None):
    df = pd.read_csv(path)
    if upper_cols:
        df.columns = df.columns.str.strip().str.upper()
    if parse_dates:
        for c in parse_dates:
            if c in df.columns:
                df[c] = pd.to_datetime(df[c], errors="coerce")
    return df

def attach_nearest_visit(pet_df, registry_df, pet_date_col, reg_date_col="EXAMDATE",
                         keep_reg_cols=("VISCODE","VISCODE2", "EXAMDATE"),
                         direction="nearest"):
    """
    Match each PET row to the nearest registry visit for the SAME RID by date.
    Requires both date cols be datetime.
    """
    pet = pet_df.copy()
    reg = registry_df.copy()

    pet = pet.dropna(subset=["RID", pet_date_col])
    reg = reg.dropna(subset=["RID", reg_date_col])

    pet = pet.sort_values(["RID", pet_date_col])
    reg = reg.sort_values(["RID", reg_date_col])

    out = pd.merge_asof(
        pet,
        reg[["RID", reg_date_col, *keep_reg_cols]],
        left_on=pet_date_col,
        right_on=reg_date_col,
        by="RID",
        direction=direction,
        suffixes=("", "_REG")
    )
    return out

def filter_pet_qc(pet_df, qc_col="QC_FLAG", qc_pass=2, tracer_col="TRACER", exclude_tracers=None):
    df = pet_df.copy()
    if qc_col in df.columns:
        df = df[df[qc_col] == qc_pass]
    if exclude_tracers and tracer_col in df.columns:
        df = df[~df[tracer_col].isin(exclude_tracers)]
    return df

def add_apoe4_status(apoe_df, genotype_col="GENOTYPE"):
    apoe = apoe_df.copy()
    apoe = apoe.drop_duplicates(subset=["RID"], keep="last")

    mapping = {"2/2":0, "2/3":0, "3/3":0, "2/4":1, "3/4":1, "4/4":2}
    apoe["APOE4"] = apoe[genotype_col].map(mapping).astype("float64")
    return apoe[["RID","APOE4"]]

# -----------------------
# 1) Load
# -----------------------
taumeta   = load_csv("YOUR/PATH/taumeta.csv", parse_dates=["SCANDATE"])
tau_qc    = load_csv("YOUR/PATH/tau_qc.csv", parse_dates=["SCANDATE"])
registry  = load_csv("YOUR/PATH/registry.csv", parse_dates=["EXAMDATE"])
tau_pet   = load_csv("YOUR/PATH/taupet.csv", parse_dates=["SCANDATE"])

adni      = load_csv("YOUR/PATH/ADNIMERGE.csv", parse_dates=["EXAMDATE","EXAMDATE_BL"])
apoe_raw  = load_csv("YOUR/PATH/apoe4.csv")

amy_meta  = load_csv("YOUR/PATH/amymeta.csv", parse_dates=["SCANDATE"])
amy_qc    = load_csv("YOUR/PATH/amyqc.csv", parse_dates=["SCANDATE"])
amy_pet   = load_csv("YOUR/PATH/amy_data.csv", parse_dates=["SCANDATE"])

fs_region = load_csv("YOUR/PATH/fs_region.csv", upper_cols=False)

# -----------------------
# 2) TAU metadata + QC (standardize key columns)
# -----------------------
# Ensure VISCODE column name exists consistently
if "VISCODE" not in tau_qc.columns and "VISCODE" in tau_qc.columns:
    pass  # (kept for clarity)

# Keep only what you need from QC
tau_qc_keep = tau_qc[["RID","VISCODE","VISCODE2","SCANDATE"]].copy()

tau = taumeta.merge(
    tau_qc_keep,
    on=["RID","VISCODE"],
    how="left",
    suffixes=("_META","_QC")
)

# Prefer taumeta scandate if present
if "SCANDATE_META" in tau.columns and "SCANDATE_QC" in tau.columns:
    tau["SCANDATE"] = tau["SCANDATE_META"].fillna(tau["SCANDATE_QC"])
    tau.drop(columns=["SCANDATE_META","SCANDATE_QC"], inplace=True, errors="ignore")

tau["SCANDATE"] = pd.to_datetime(tau["SCANDATE"], errors="coerce")

# -----------------------
# 3) Attach nearest registry visit to TAU (CORRECT: by RID)
# -----------------------
tau_with_visit = attach_nearest_visit(tau, registry, pet_date_col="SCANDATE")

# -----------------------
# 4) Merge TAU PET data onto visit-annotated TAU
# -----------------------
tau_pet = filter_pet_qc(tau_pet, qc_col="QC_FLAG", qc_pass=2, exclude_tracers=["MK6240"])
tau_merged = tau_pet.merge(
    tau_with_visit,
    on=["RID","VISCODE"],
    how="left",
    suffixes=("", "_TAUINFO")
)
tau_merged = tau_merged.loc[:, ~tau_merged.columns.duplicated()].copy()

# -----------------------
# 5) Attach nearest ADNIMERGE row to TAU (by RID)
# -----------------------
# Select columns you want (keep RID, EXAMDATE, VISCODE at minimum)
adni_keep_cols = ["RID","VISCODE","EXAMDATE","EXAMDATE_BL","AGE","DX"]  # + your extras
adni2 = adni[[c for c in adni_keep_cols if c in adni.columns]].copy()

tau_adni = pd.merge_asof(
    tau_merged.sort_values(["RID","SCANDATE"]),
    adni2.sort_values(["RID","EXAMDATE"]),
    left_on="SCANDATE",
    right_on="EXAMDATE",
    by="RID",
    direction="nearest",
    suffixes=("", "_ADNI")
)

# -----------------------
# 6) APOE4 status + Age at scan
# -----------------------
apoe = add_apoe4_status(apoe_raw, genotype_col="GENOTYPE")
tau_adni = tau_adni.merge(apoe, on="RID", how="left")

# age at scan requires AGE (at EXAMDATE) + offset from EXAMDATE_BL
if {"AGE","EXAMDATE_BL","SCANDATE"}.issubset(tau_adni.columns):
    tau_adni["AGE_AT_SCAN"] = tau_adni["AGE"] + (tau_adni["SCANDATE"] - tau_adni["EXAMDATE_BL"]).dt.days / 365.25

tau_adni = tau_adni.dropna(subset=["APOE4","DX"])

# -----------------------
# 7) Amyloid meta + QC + registry visit
# -----------------------
amy = amy_meta.merge(
    amy_qc[["RID","SCANDATE","VISCODE","VISCODE2"]].copy(),
    on="RID",
    how="left",
    suffixes=("_META","_QC")
)

# Prefer meta date if present
if "SCANDATE_META" in amy.columns and "SCANDATE_QC" in amy.columns:
    amy["SCANDATE"] = amy["SCANDATE_META"].fillna(amy["SCANDATE_QC"])
    amy.drop(columns=["SCANDATE_META","SCANDATE_QC"], inplace=True, errors="ignore")

amy["SCANDATE"] = pd.to_datetime(amy["SCANDATE"], errors="coerce")
amy = amy.dropna(subset=["RID","SCANDATE"]).drop_duplicates(subset=["RID","SCANDATE"])

amy_with_visit = attach_nearest_visit(amy, registry, pet_date_col="SCANDATE")

# -----------------------
# 8) Merge amyloid PET + QC
# -----------------------
amy_pet = filter_pet_qc(amy_pet, qc_col="QC_FLAG", qc_pass=2)
amy_pet = amy_pet.dropna(subset=["RID","VISCODE","SCANDATE"])

amy_merged = amy_with_visit.merge(
    amy_pet,
    on=["RID","VISCODE"],
    how="inner",
    suffixes=("", "_AMYPET")
)

# -----------------------
# 9) Merge tau + amyloid by nearest scan date (by RID) and threshold
# -----------------------
tau_adni = tau_adni.rename(columns={"SCANDATE":"SCANDATE_TAU"})
amy_merged = amy_merged.rename(columns={"SCANDATE":"SCANDATE_AMY"})

tau_amy = pd.merge_asof(
    tau_adni.sort_values(["RID","SCANDATE_TAU"]),
    amy_merged.sort_values(["RID","SCANDATE_AMY"]),
    left_on="SCANDATE_TAU",
    right_on="SCANDATE_AMY",
    by="RID",
    direction="nearest",
    suffixes=("_TAU","_AMY")
)

tau_amy["AMY_TAU_SCANDIFF_DAYS"] = (tau_amy["SCANDATE_AMY"] - tau_amy["SCANDATE_TAU"]).dt.days
tau_amy = tau_amy[tau_amy["AMY_TAU_SCANDIFF_DAYS"].abs() < 365]

# -----------------------
# 10) FS region matching (on the FINAL dataset)
# -----------------------
fs_region["FS_LABEL"] = fs_region["FS_LABEL"].str.replace("-","_", regex=False).str.upper()
fs_region["FS_LABEL"] = fs_region["FS_LABEL"] + "_SUVR"

suvr_cols = [c for c in fs_region["FS_LABEL"].tolist() if c in tau_amy.columns]

# drop subcortical indices if your fs_region file matches that ordering
subcortical_idx = [34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,83]
ordered = np.delete(np.array(fs_region["FS_LABEL"].tolist(), dtype=object), subcortical_idx).tolist()
suvr_cols_ordered = [c for c in ordered if c in tau_amy.columns]

non_suvr_cols = [c for c in tau_amy.columns if not c.endswith("_SUVR")]
tau_amy_cortical = tau_amy[non_suvr_cols + suvr_cols_ordered].copy()
tau_amy_cortical.to_csv("final_matched_dataset.csv", index = False)

print("Final matched dataset:", tau_amy_cortical.shape)
