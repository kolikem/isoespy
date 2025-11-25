""" 目的
    1. hmmscan の結果から有意なドメインを回収
    2. transcript-based annotation のフォーマット（GTF化の前段階）で出力
"""

import pandas as pd
import argparse
import sys
import csv

# #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
# #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
# 7tm_2                PF00002.29   250 ENST00000283303.3    -            695   2.9e-33  115.5  23.6   1   1   6.9e-37   5.5e-33  114.7  23.6     3   250   400   648   398   648 0.86 7 transmembrane receptor (Secretin family)

# |
# |
# |
# V

# #transcript_id	from	to	nt/aa		feature_name		source		attribute_field
# ADGRF4_novel_tx_1	252     500	amino_acid	protein_domain (引数)	Pfam (引数)	domain_id "PF22691.2"; domain_name "Thiolase_C_1"; evalue "4.3e-36";


# hmmscan の domtbl をパース
def parse_hmmscan_result(file, coord):
    if coord == "ali":	# ali
        C = 17
    else:		# env
        C = 19

    df = []
    with open(file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            df.append([
                cols[3],			# 転写産物名
                cols[0],			# PFAM shortname
                cols[1],			# PFAM accession
                float(cols[12]),		# i-Evalue
                int(cols[2]),			# target length
                int(cols[15]), int(cols[16]),	# hmm from, hmm to
                int(cols[C]), int(cols[C+1]),	# tx from, tx to
                float(cols[21]),		# The mean posterior probability of aligned residues in the MEA alignment
                int(cols[5]),			# CDS total length
            ])
    return pd.DataFrame(df, columns=["Transcript", "Shortname", "Pfam_accession", "Evalue", "Domain_len", "hmm_from", "hmm_to", "tx_from", "tx_to", "Accuracy", "CDS_total_len"])


def remove_overlapping_domains(df):
    result_list = []
    for transcript, group in df.groupby("Transcript"):
        group = group.sort_values("Evalue").reset_index(drop=True)
        kept = []

        for i, row in group.iterrows():
            # 現在の領域
            from_i, to_i = row["tx_from"], row["tx_to"]
            overlap = False

            for kept_row in kept:
                from_k, to_k = kept_row["tx_from"], kept_row["tx_to"]
                # 重なり判定（両端含む）
                if not (to_i < from_k or from_i > to_k):
                    overlap = True
                    break

            if not overlap:
                kept.append(row)

        result_list.extend(kept)

    return pd.DataFrame(result_list)


def save_domain_annotation_tsv(filtered_df, output_path, feature, source):
    """
    指定された filtered_df からドメイン注釈ファイルを作成して保存する。

    Parameters:
        filtered_df (pd.DataFrame): ドメイン予測結果のフィルタ済みデータフレーム
        output_path (str): 出力ファイルパス（TSV）
        feature (str): feature 列に使う固定文字列
        source (str): source 列に使う固定文字列
    """
    # attributes列を作成
    def format_attributes(row):
        return (
            'domain_id ' + '"'+str(row["Pfam_accession"]).replace('"', '')+'"' + '; '
            'domain_name ' + '"'+str(row["Shortname"]).replace('"', '')+'"' + '; '
            'evalue ' + '"'+str(row["Evalue"]).replace('"', '')+'"' + ';'
        )

    attributes = filtered_df.apply(format_attributes, axis=1)

    output_df = pd.DataFrame({
        "#Transcript": filtered_df["Transcript"],
        "tx_from": filtered_df["tx_from"],
        "tx_to": filtered_df["tx_to"],
        "amino_acid": "amino_acid",
        "feature": feature,
        "source": source,
        "attributes": attributes
    })

    # 保存
    output_df.to_csv(output_path, sep="\t", index=False, quoting=csv.QUOTE_NONE)


def filter_hmmscan_result(df, evalue_threshold, accuracy_threshold, coverage_threshold, output_path, feature, source, overlap_remover):
    """
    HMMER hmmscanの出力データをフィルタリングし、条件を満たすドメインのみを保存する。

    Parameters:
        evalue_threshold (float): E-valueのカットオフ（デフォルト: 1e-5）
        accuracy_threshold (float): Accuracyのカットオフ（デフォルト: 0.7）
        coverage_threshold (float): アライメントカバレッジのカットオフ（デフォルト: 0.5）
    """

    # アライメント割合を計算
    df["alignment_ratio"] = (df["hmm_to"] - df["hmm_from"] + 1) / df["Domain_len"]

    # (1) Evalue < evalue_threshold
    filtered_df = df[df["Evalue"] < evalue_threshold].copy()

    # (2) アライメント割合 > alignment_ratio_threshold
    filtered_df = filtered_df[filtered_df["alignment_ratio"] > coverage_threshold]

    # (3) Accuracy > accuracy_threshold
    filtered_df = filtered_df[filtered_df["Accuracy"] > accuracy_threshold]

    # (4) (overlap_remover == True の場合のみ) Transcriptごとに予測領域の重なるドメインはEvalueが最小のものに絞る
    if overlap_remover:
        filtered_df = remove_overlapping_domains(filtered_df)

    # 結果を保存
    output_path = output_path+".gtfprep"
    save_domain_annotation_tsv(filtered_df, output_path, feature, source)


def str2bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() == "true":
        return True
    elif value.lower() == "false":
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected: 'true' or 'false'")


def main():
    parser = argparse.ArgumentParser(description="Filter hmmscan results")
    parser.add_argument("domtbl", help="Path to domtbl file")
    parser.add_argument("feature", help="feature name column")
    parser.add_argument("source", help="source name column")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold (default: 1e-5)")
    parser.add_argument("--coverage", type=float, default=0.5, help="Minimum alignment coverage ratio of the query to the full length of  domain model (default: 0.5)")
    parser.add_argument("--acc", type=float, default=0.7, help="Accuracy threshold (default: 0.7)")
    parser.add_argument("--filt_overlap", type=str2bool, default=False, help="true or false. true if overlapped domains are not allowed. (default: False)")

    args = parser.parse_args()

    # 必須引数
    domtbl = args.domtbl
    feature = args.feature
    source = args.source

    # オプション
    evalue_threshold = args.evalue
    coverage_threshold = args.coverage
    acc_threshold = args.acc
    coord = "ali"
    overlap_remover = args.filt_overlap

    # メイン処理
    df = parse_hmmscan_result(domtbl, coord)
    filter_hmmscan_result(df, evalue_threshold, acc_threshold, coverage_threshold, domtbl, feature, source, overlap_remover)

if __name__ == "__main__":
    main()

