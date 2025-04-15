""" 目的
    1. SignalP6 の結果を回収
    2. GTFprepフォーマットで出力
"""

import pandas as pd
import argparse
import sys
import csv

### output.gff3
#       ## gff-version 3
#       ENST00000261405.10      SignalP-6.0     signal_peptide  1       22      0.99975723      .       .       .
#       ENST00000216336.3       SignalP-6.0     signal_peptide  1       18      0.99973375      .       .       .

# |
# |
# |
# V

# #transcript_id	from	to	nt/aa		feature_name		source		attribute_field
# ENST00000261405.10	1	22	amino_acid	signal_peptide		SignalP-6.0	likelihood "0.99975723";


def parse_result(file):
    df = []
    with open(file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            df.append([
                cols[0],	# 転写産物名
                int(cols[3]),	# from
                int(cols[4]),	# to
                float(cols[5]),	# likelihood
            ])
    return pd.DataFrame(df, columns=["Transcript", "from", "to", "likelihood"])


def filter_overlaped(df):
    df_filtered = df.loc[df.groupby('Transcript')['likelihood'].idxmax()]
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
            'likelihood ' + '"'+str(row["likelihood"]).replace('"', '')+'"' + '; '
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

def main():
    parser = argparse.ArgumentParser(description="Make gtfprep format from signalp results")
    parser.add_argument("signalp_output_gff3", help="Path to output.gff3")
    parser.add_argument("feature", help="feature name column")
    parser.add_argument("source", help="source name column")

    args = parser.parse_args()

    # 必須引数
    gff3 = args.signalp_output_gff3
    feature = args.feature
    source = args.source

    # メイン処理
    output_path = gff3 + ".gtfprep"
    df = parse_result(gff3)
    filtered_df = filter_overlaped(df)
    save_gtfprep(filtered_df, output_path, feature, source)

if __name__ == "__main__":
    main()

