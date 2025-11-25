import shutil
import subprocess
import sys
import argparse
import os
import csv
import re
import pandas as pd
from collections import OrderedDict

def run_r_script(
    r_script_path,
    input1,
    input2,
    input3,
    output1,
    output2
):
    """
    Run an R script via subprocess with 3 input files and 2 output files.

    Parameters
    ----------
    r_script_path : str
        Path to the R script (e.g., "script_A.R")
    input1, input2, input3 : str
        Paths to input files
    output1, output2 : str
        Paths to output files
    """

    # check that Rscript is available
    if not shutil.which("Rscript"):
        raise EnvironmentError("Rscript not found in PATH. Please ensure R is installed and accessible.")

    # construct command
    cmd = [
        "Rscript",
        r_script_path,
        input1,
        input2,
        input3,
        output1,
        output2
    ]

    print("Running command:")
    print(" ".join(cmd))

    try:
        subprocess.run(cmd, check=True)
        print("R script executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: R script failed with exit code {e.returncode}")
        sys.exit(1)


def parse_attributes(attr_field):
    """
      - key "value"; key "value";
      - key=value;key=value;
      - key value; のような形式
    """
    attr = {}
    if attr_field is None:
        return attr
    s = attr_field.strip().strip(";")
    if not s:
        return attr

    # ; で分割（連続セミコロンや末尾空要素を無視）
    parts = re.split(r';\s*', s)
    for part in parts:
        part = part.strip()
        if not part:
            continue
        key = None
        val = None

        if '=' in part:
            key, val = part.split('=', 1)
        else:
            # key "value" あるいは key value を拾う
            m = re.match(r'(\S+)\s+"?(.*?)"?$', part)
            if m:
                key, val = m.group(1), m.group(2)

        if key is not None:
            key = key.strip()
            if val is not None:
                val = val.strip().strip('"')
            attr[key] = val
    return attr

def build_tx_to_gene_map(gtf_tsv_path, gene_key, tx_key):
    """
    第9列 attributes から gene_key と tx_key を取り出し、
    転写産物ID -> 遺伝子ID の辞書を作る。
    同一転写産物IDが複数行に出ても最初の対応を採用（2回目以降は無視）。
    """
    tx2gene = OrderedDict()
    line_no = 0
    with open(gtf_tsv_path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            line_no += 1
            if not row or len(row) < 9:
                # 不完全行はスキップ
                continue
            attributes_field = row[8]
            attrs = parse_attributes(attributes_field)

            tx = attrs.get(tx_key)
            gene = attrs.get(gene_key)

            if tx is None or gene is None:
                # 必要キーが無い行はスキップ（構造上、exon行など全てに両方あるとは限らないため）
                continue

            if tx not in tx2gene:
                tx2gene[tx] = gene
            else:
                # 2回目以降は飛ばす（ただし gene が不一致なら警告）
                if tx2gene[tx] != gene:
                    print(
                        f"[WARN] Transcript {tx!r} already mapped to {tx2gene[tx]!r}, "
                        f"but found another gene {gene!r} at line {line_no}. Keeping the first mapping.",
                        file=sys.stderr
                    )
    if not tx2gene:
        print("[ERROR] No (transcript_id, gene_id) pairs were extracted from attributes. "
              f"Check that 9th column contains attributes and keys '{tx_key}', '{gene_key}'.",
              file=sys.stderr)
        sys.exit(1)
    return tx2gene

def read_counts(counts_path):
    """
    counts TSV を読み込み：
      - 1行目: ヘッダー（col1=tx、col2+=サンプル名）
      - 2行目以降: tx, sample1_count, sample2_count, ...
    戻り値:
      header_samples: [sample1, sample2, ...]
      rows: [(tx, [c1, c2, ...]), ...]  ※入力の順序維持
    """
    rows = []
    with open(counts_path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            print("[ERROR] Counts file is empty.", file=sys.stderr)
            sys.exit(1)
        if len(header) < 2:
            print("[ERROR] Counts header must have at least two columns: tx and >=1 sample.", file=sys.stderr)
            sys.exit(1)

        # 1列目は tx、2列目以降はサンプル名
        sample_headers = header[1:]

        for i, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) < 2:
                print(f"[ERROR] Counts row {i} has fewer than 2 columns.", file=sys.stderr)
                sys.exit(1)
            tx = row[0]
            counts = row[1:]
            rows.append((tx, counts))

    if not rows:
        print("[ERROR] Counts file has no data rows.", file=sys.stderr)
        sys.exit(1)

    return sample_headers, rows


def ensure_sets_match(tx2gene, counts_rows):
    """
    GTF は発現量データをグループ化するために使うだけなので、
    GTF にしかない transcript は許容する。
    ただし、counts 側にある transcript が GTF に無いのはエラーにする。
    """
    gtf_txs    = set(tx2gene.keys())               # GTFに書いてあるtranscript
    counts_txs = set(tx for tx, _ in counts_rows)  # countsに登場するtranscript

    # counts にあるのに GTF に無いものだけチェック
    missing_in_gtf = sorted(counts_txs - gtf_txs)

    if missing_in_gtf:
        msg_lines = ["[ERROR] Some transcripts in counts were not found in GTF (these cannot be grouped by gene)."]
        msg_lines.append(f"  - Missing in GTF: {len(missing_in_gtf)} examples: {missing_in_gtf[:10]}")
        msg_lines.append("Please make sure every transcript in the expression data exists in the GTF.")
        print("\n".join(msg_lines), file=sys.stderr)
        sys.exit(1)

    # 逆に「GTFにだけあるもの」は何もしないでOK


def write_output(output_path, sample_headers, counts_rows, tx2gene):
    """
    DRIMSeq 形式のTSV:
      gene_id  feature_id  sample1  sample2  ...
      (fileBの行順を維持)
    """
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        header = ["gene_id", "feature_id"] + sample_headers
        writer.writerow(header)

        for tx, counts in counts_rows:
            gene = tx2gene.get(tx)
            if gene is None:
                # 通常ここには来ない（事前一致チェック済み）
                print(f"[WARN] tx {tx!r} not found in mapping; skipping.", file=sys.stderr)
                continue
            writer.writerow([gene, tx] + counts)


def merge_results(pval_file, prop_file, sample_file, output_file):
    NORM = "dot2dash"
    # --- load ---
    df_stage = pd.read_csv(pval_file, sep="\t")
    df_prop  = pd.read_csv(prop_file, sep="\t")
    df_smp   = pd.read_csv(sample_file, sep="\t")

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
        if NORM == "dot2dash":
            return s.replace(".", "-")
        elif NORM == "dash2dot":
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
    out.to_csv(output_file, sep="\t", index=False)
    print(f"Done: {output_file}")


def main():
    # Arguments
    parser = argparse.ArgumentParser(description="Run an R script with 3 input and 2 output files.")
    parser.add_argument("--r_script", required=True, help="Path to the R script (e.g., script_A.R)")
    parser.add_argument("--count", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output_pval", required=True)
    parser.add_argument("--output_prop", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--gene_key", default="gene_symbol", help="Attribute key for gene id (default: gene_symbol)")
    parser.add_argument("--tx_key", default="transcript_id", help="Attribute key for transcript id (default: transcript_id)")
    parser.add_argument("--OUT", default="result.tsv", help="Output TSV")

    args = parser.parse_args()

    # Prepare DRIMSeq input files
    tx2gene = build_tx_to_gene_map(args.gtf, args.gene_key, args.tx_key)
    sample_headers, counts_rows = read_counts(args.count)
    ensure_sets_match(tx2gene, counts_rows)
    drimseq_count = os.path.dirname(args.sample)+"/"+os.path.basename(args.count)+".DRIMSeq.tsv"
    write_output(drimseq_count, sample_headers, counts_rows, tx2gene)
    print(f"[OK] Wrote DRIMSeq counts table: {drimseq_count}", file=sys.stderr)

    # DRIMSeq + stageR
    run_r_script(
        args.r_script,
        drimseq_count,
        args.sample,
        args.config,
        args.output_pval,
        args.output_prop
    )

    # Merge results
    merge_results(args.output_pval, args.output_prop, args.sample, args.OUT)


if __name__ == "__main__":
    main()


