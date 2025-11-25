""" 転写産物上のアミノ酸範囲 ==> 転写産物上の塩基配列範囲 ==> ゲノムエクソン上の範囲
　　に変換するモジュール
"""

# ================= Input =================
# arg1 - GTF (CDS feature 必須)
# arg2 - 転写産物名
# arg3 - アミノ酸範囲 (1-based, inclusive)
# arg4 - cds全長

# ================= Output =================
# ゲノムにおける範囲 (list)　e.g., [(start1, end1), (start2, end2), (start3, end3)]    **start昇順


def aanucConverter(aafrom, aato, frame):
    """ アミノ酸配列上範囲 ==> CDS座標上範囲 """
    nucfrom = 3*aafrom - 2 + frame
    nucto = 3*aato + frame
    return (nucfrom, nucto)

def map_feature_to_genome(coord_tx, coord_gen, strand_gen):
    """ 塩基配列上範囲 ==> ゲノム上範囲 """

    start_tx = coord_tx["feature_start"]	# 転写産物塩基配列上のfeature 開始点 CDS (aa), exon (nuc)
    end_tx = coord_tx["feature_end"]		# 転写産物塩基配列上のfeature 終了点
    totallen_tx = coord_tx["total_len"]		# CDS全長 (aa) / exon全長 (nuc)

    if strand_gen == "+":   # ゲノム上でアイソフォームが+鎖にコード
        start_tx, end_tx = start_tx, end_tx
    else:                   # ゲノム上でアイソフォームが-鎖にコード
        start_tx, end_tx = totallen_tx - end_tx + 1, totallen_tx - start_tx + 1

    exon_table = {"tx": [], "gen": []}		# tx: 転写産物上位置　　gen: ゲノム上位置
    feature_table = {"gen": []}			# ゲノム上 feature 位置

    # exon 位置対応表 tx --- genome を作成
    exon_table["gen"] = coord_gen
    cds_pos = 0
    for exon_i in exon_table["gen"]:
        exon_len = exon_i[1] - exon_i[0] + 1
        exon_table["tx"].append((cds_pos+1, cds_pos+exon_len))
        cds_pos += exon_len

    #print("exon_table", exon_table)
    if totallen_tx != exon_table["tx"][-1][1]:
        print("この転写産物は全長が不一致です!!!")

    #print(exon_table)

    # feature 座標をゲノムにマップする
    for i in range(len(exon_table["tx"])):
        if exon_table["tx"][i][0] <= start_tx <= exon_table["tx"][i][1]:
            start_exon = i + 1									# start_exon: feature が開始するエクソン番号 (参照ゲノム+鎖上流から数えて)
            start_pos = exon_table["gen"][i][0] + (start_tx - exon_table["tx"][i][0])		# start_pos: ゲノム上のfeature 開始点
        if exon_table["tx"][i][0] <= end_tx <= exon_table["tx"][i][1]:
            end_exon = i + 1									# end_exon: feature が終了するエクソン番号 (参照ゲノム+鎖上流から数えて)
            end_pos = exon_table["gen"][i][0] + (end_tx - exon_table["tx"][i][0])		# end_pos: ゲノム上のfeature 終了点

    if start_exon == end_exon:   # feature が１つのエクソン内で完結
        feature_table["gen"].append((start_pos, end_pos))
    else:                        # feature が複数のエクソンにまたがる
        for i in range(start_exon-1, end_exon):
            if i == start_exon-1:   # CDSを含む始めのエクソン
                feature_table["gen"].append((start_pos, exon_table["gen"][i][1]))
            elif i == end_exon-1:   # CDSを含む最後のエクソン
                feature_table["gen"].append((exon_table["gen"][i][0], end_pos))
            else:                   # CDSを含む内部のエクソン
                feature_table["gen"].append(exon_table["gen"][i])

    #print("feature_table")
    #print(feature_table)
    return feature_table["gen"]   # タプルリストを返す


def backmap_aa(coord_gen, tx, aafrom, aato, cds_total_len, strand_gen, frame_d):
    """ coord_gen : "CDS" のゲノム座標上での位置情報 list形式
        tx        : 転写産物名
        aafrom    : アミノ酸配列上におけるfeature (ドメインなど) の開始　  1-based, inclusive
        aato      : アミノ酸配列上におけるfeature (ドメインなど) の終了　　1-based, inclusive
    """
    if tx in frame_d:
        frame = frame_d[tx]
    else:
        # フレームデータが無い時は frame = 0 と仮定
        frame = 0

    # アミノ酸配列上位置 ==> 塩基配列上配列位置 (CDS座標上, 転写産物座標上ではない)
    nucfrom, nucto = aanucConverter(aafrom, aato, frame)

    # CDS座標上位置 ==> ゲノム上位置
    coord_cds = {"feature_start": nucfrom, "feature_end": nucto, "total_len": cds_total_len}
    backmap_coords = map_feature_to_genome(coord_cds, coord_gen, strand_gen)
    return backmap_coords


def backmap_nuc(coord_gen, tx, nucfrom, nucto, exon_total_len, strand_gen):
    """ coord_gen : "exon" のゲノム座標上での位置情報 list形式
        tx        : 転写産物名
        nucfrom   : 塩基配列上におけるfeature (CDSなど) の開始　  1-based, inclusive
        nucto     : 塩基配列上におけるfeature (CDSなど) の終了　　1-based, inclusive
    """
    # 転写産物上位置 ==> ゲノム上位置
    coord_exon = {"feature_start": nucfrom, "feature_end": nucto, "total_len": exon_total_len}
    backmap_coords = map_feature_to_genome(coord_exon, coord_gen, strand_gen)
    return backmap_coords


if __name__ == "__main__":
    print("Call from another script.")

