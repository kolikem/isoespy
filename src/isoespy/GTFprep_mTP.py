""" 目的
    1. TargetP の結果を回収
    2. GTFprepフォーマットで出力
"""

import re
import pandas as pd
import argparse
import sys
import csv


### TargetP2 (mTPを含む)
#       # TargetP-2.0   Organism: Non-Plant     Timestamp: 20250314183052
#       # ID    		Prediction      noTP    	SP      	mTP     	CS Position
#       ENST00000465205.2       mTP     	0.296826        0.000283        0.702891        CS pos: 25-26. VMG-KQ. Pr: 0.1151
#       ENST00000261405.10      SP      	0.000257        0.999722        0.000021        CS pos: 22-23. TLC-AE. Pr: 0.8440


def get_cleavage_site(text):
    # 正規表現でCS posの2つの数字を取り出す
    match = re.search(r'CS pos:\s*(\d+)-(\d+)', text)
    if match:
        cs_start = int(match.group(1))
        cs_end = int(match.group(2))
    return cs_start


def get_csprob(text):
    # 正規表現でCS posの確率を取り出す
    match = re.search(r'Pr:\s*([0-9]*\.[0-9]+)', text)
    if match:
        cs_prob = float(match.group(1))
    return cs_prob


def parse_result(file):
    df = []
    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue

            if line.startswith("<query>"):
                continue

            if not line:
                continue

            cols = line.split('\t')
            if len(cols) != 6:
                continue

            if cols[1] != "mTP":
            # mTP prediction のみ処理
                continue

            if "?" in cols[5]:
                continue
            # Cleavage site
            cleavage_site = get_cleavage_site(cols[5])
            cleavage_prob = get_csprob(cols[5])

            df.append([
                cols[0],		# 転写産物名
                1,			# from
                cleavage_site,		# to
                float(cols[4]),		# mTP prob.
                cleavage_prob,		# cleavage site prob.
            ])
    return pd.DataFrame(df, columns=["Transcript", "from", "to", "mTP_prob", "cleavage_prob"])


def filter_threshold(df, CUTOFF, SINGLE):
    # CUTOFF 未満の結果をフィルター（デフォルト：0.9）
    df_filtered = df[df["mTP_prob"] >= CUTOFF]

    # 転写産物につき1つの結果だけを残す（デフォルト：しない）
    df_filtered = df_filtered.loc[df_filtered.groupby('Transcript')['cleavage_prob'].idxmax()]
    return df_filtered


def save_gtfprep(filtered_df, output_path, feature, source):
    """
    Parameters:
        filtered_df (pd.DataFrame): フィルタ済みデータフレーム
        output_path (str): 出力ファイルパス（TSV）
        feature (str): feature 列に使う固定文字列
        source (str): source 列に使う固定文字列
    """
    # attributes列を作成
    def format_attributes(row):
        return (
            'mTP ' + '"'+str(row["mTP_prob"]).replace('"', '')+'"' + '; '
            'CS_prob ' + '"'+str(row["cleavage_prob"]).replace('"', '')+'"' + '; '
        )

    attributes = filtered_df.apply(format_attributes, axis=1)

    output_df = pd.DataFrame({
        "#Transcript": filtered_df["Transcript"],
        "tx_from": filtered_df["from"],
        "tx_to": filtered_df["to"],
        "amino_acid": "amino_acid",
        "feature": feature,
        "source": source,
        "attributes": attributes
    })

    # 保存
    output_df.to_csv(output_path, sep="\t", index=False, quoting=csv.QUOTE_NONE)
    print(">>> Output file: ", output_path)


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
    parser = argparse.ArgumentParser(description="Make gtfprep format from TargetP results")
    parser.add_argument("output", help="Path to TargetP output file")
    parser.add_argument("feature", help="feature name column")
    parser.add_argument("source", help="source name column")
    parser.add_argument("--cutoff", type=float, default=0.9, help="Cutoff of mTP prob (default: 0.9)")
    parser.add_argument("--single", type=str2bool, default=False, help="true or false. True if you keep only one detected result for each transcript. (default: False)")

    args = parser.parse_args()

    # 引数
    nls_output = args.output
    feature = args.feature
    source = args.source
    cutoff = args.cutoff
    single = args.single

    # メイン処理
    output_path = nls_output + ".gtfprep"
    df = parse_result(nls_output)
    filtered_df = filter_threshold(df, cutoff, single)
    save_gtfprep(filtered_df, output_path, feature, source)

if __name__ == "__main__":
    main()

