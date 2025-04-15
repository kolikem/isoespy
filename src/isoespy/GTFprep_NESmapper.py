""" 目的
    1. NESmapper の結果を回収
    2. GTFprepフォーマットで出力
"""

import pandas as pd
import argparse
import sys
import csv

### NESmapper
#	<query>			<pos>	<nes>			<score>		<nes-class>
#	AADAT_novel_tx_1	216	IYELARKYDFLIIE		18		class_1a
#	ABCC2_novel_tx_1	566	RPELDLVLRGITCD		12.6		class_1a
#	ABCE1_novel_tx_1	27	KSCPVVRMGKLCIE		12.6		class_1b


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
            if len(cols) != 5:
                continue

            df.append([
                cols[0],			# 転写産物名
                int(cols[1]),			# from
                int(cols[1])+len(cols[2])-1,	# to
                float(cols[3]),			# score
                cols[4],			# nes-class
            ])
    return pd.DataFrame(df, columns=["Transcript", "from", "to", "score", "nes-class"])


def filter_threshold(df, CUTOFF, SINGLE):
    # CUTOFF 未満の結果をフィルター（デフォルト：12）
    df_filtered = df[df["score"] >= CUTOFF]

    # 転写産物につき1つの結果だけを残す（デフォルト：しない）
    df_filtered = df_filtered.loc[df_filtered.groupby('Transcript')['score'].idxmax()]
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
            'score ' + '"'+str(row["score"]).replace('"', '')+'"' + '; '
            'new-class ' + '"'+str(row["nes-class"]).replace('"', '')+'"' + '; '
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
    parser = argparse.ArgumentParser(description="Make gtfprep format from NESmapper results")
    parser.add_argument("output", help="Path to NESmapper output file")
    parser.add_argument("feature", help="feature name column")
    parser.add_argument("source", help="source name column")
    parser.add_argument("--cutoff", type=float, default=12, help="Cutoff of score (default: 12)")
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

