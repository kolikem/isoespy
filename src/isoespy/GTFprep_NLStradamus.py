""" 目的
    1. NLStradamus の結果を回収
    2. GTFprepフォーマットで出力
"""

import pandas as pd
import argparse
import sys
import csv


### NLStradamus
#	ENST00000358030.6	posterior	0.920	176	199	RPRKKGMSLKSDAVLSKSKRRKKP
#	Finished analyzing ENST00000358030.6. Found 1 sites.
#
#	ENST00000345306.10	posterior	0.920	341	364	RPRKKGMSLKSDAVLSKSKRRKKP
#	Finished analyzing ENST00000345306.10. Found 1 sites.
#
#	===================================================
#	Analyzed 26 proteins.
#	2 sites were found using the posterior probability threshold.
#	Input file : trial2/test.TransDecoder.gtf.translate.faa.
#	Threshold used : 0.6.
#	===================================================


def parse_result(file):
    df = []
    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue

            if line.startswith("Finished"):
                continue

            if not line:
                continue

            if  line.startswith("==="):
                break

            cols = line.split('\t')

            df.append([
                cols[0],	# 転写産物名
                int(cols[3]),	# from
                int(cols[4]),	# to
                float(cols[2]),	# posterior
            ])
    return pd.DataFrame(df, columns=["Transcript", "from", "to", "posterior"])


def filter_threshold(df, CUTOFF, SINGLE):
    # CUTOFF 未満の結果をフィルター（デフォルト：0.9）
    df_filtered = df[df["posterior"] >= CUTOFF]

    # 転写産物につき1つの結果だけを残す（デフォルト：しない）
    df_filtered = df_filtered.loc[df_filtered.groupby('Transcript')['posterior'].idxmax()]
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
            'posterior ' + '"'+str(row["posterior"]).replace('"', '')+'"' + '; '
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
    parser = argparse.ArgumentParser(description="Make gtfprep format from NLStradamus results")
    parser.add_argument("output", help="Path to NLStradamus output file")
    parser.add_argument("feature", help="feature name column")
    parser.add_argument("source", help="source name column")
    parser.add_argument("--cutoff", type=float, default=0.9, help="Cutoff of posterior probability (default: 0.9)")
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

