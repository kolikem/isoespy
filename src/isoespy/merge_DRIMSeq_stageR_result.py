#!/usr/bin/env python3
import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser(
    description="Append group-wise mean transcript proportions to stageR TSV."
)
parser.add_argument("--stageR", required=True)
parser.add_argument("--proportions", required=True)
parser.add_argument("--samples", required=True)
parser.add_argument("--out", default="stageR_with_meanProportions.tsv")
parser.add_argument("--id-normalize", choices=["none", "dot2dash", "dash2dot"],
                    default="dot2dash",
                    help="Normalize sample IDs: dot2dash (RK002.N->RK002-N), "
                         "dash2dot (RK002-N->RK002.N), or none.")
args = parser.parse_args()

# --- load ---
df_stage = pd.read_csv(args.stageR, sep="\t")
df_prop  = pd.read_csv(args.proportions, sep="\t")
df_smp   = pd.read_csv(args.samples, sep="\t")

# --- rename txID->feature_id if needed ---
if "txID" in df_stage.columns:
    df_stage = df_stage.rename(columns={"txID": "feature_id"})

# --- choose id columns and sample columns ---
id_cols = [c for c in ["gene_id", "feature_id"] if c in df_prop.columns]
if not set(["gene_id","feature_id"]).issubset(df_prop.columns):
    sys.exit("ERROR: proportions TSV must contain gene_id and feature_id columns.")

sample_cols = [c for c in df_prop.columns if c not in id_cols]
if len(sample_cols) == 0:
    sys.exit("ERROR: No sample columns detected in proportions TSV.")

# --- normalize sample IDs consistently with samples.tsv ---
def normalize(s: str) -> str:
    if args.id_normalize == "dot2dash":
        return s.replace(".", "-")
    elif args.id_normalize == "dash2dot":
        return s.replace("-", ".")
    return s

# make a mapping original_sample_col -> normalized sample_id
norm_map = {col: normalize(col) for col in sample_cols}

# melt to long
df_long = df_prop.melt(
    id_vars=id_cols,
    value_vars=sample_cols,
    var_name="sample_col",
    value_name="proportion"
)
# attach normalized sample_id
df_long["sample_id"] = df_long["sample_col"].map(norm_map)

# --- join group info ---
if not set(["sample_id","group"]).issubset(df_smp.columns):
    sys.exit("ERROR: samples TSV must contain 'sample_id' and 'group' columns.")
df_long = df_long.merge(df_smp[["sample_id","group"]], on="sample_id", how="left")

# quick sanity check: unmatched samples
n_unmatched = df_long["group"].isna().sum()
if n_unmatched > 0:
    uniq_unmatched = sorted(df_long.loc[df_long["group"].isna(), "sample_id"].unique().tolist())
    print(f"[WARN] {n_unmatched} rows had no group (first few unmatched IDs: {uniq_unmatched[:10]})", file=sys.stderr)

# --- compute group-wise mean proportions ---
grp_cols = id_cols + ["group"]
mean_prop = (df_long.dropna(subset=["group"])  # drop rows without group
             .groupby(grp_cols, as_index=False)["proportion"]
             .mean())

# pivot to wide: prop_<group>
wide = (mean_prop
        .pivot_table(index=id_cols, columns="group", values="proportion")
        .reset_index())
wide = wide.rename(columns={c: f"prop_{c}" for c in wide.columns if c not in id_cols})

# --- merge to stageR (by feature_id) ---
if "feature_id" not in df_stage.columns:
    sys.exit("ERROR: stageR TSV must contain column 'feature_id' (or txID).")
out = df_stage.merge(wide, on="feature_id", how="left")

# --- write ---
out.to_csv(args.out, sep="\t", index=False)
print(f"Done: {args.out}")
