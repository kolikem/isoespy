import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as mtransforms
import copy
import pandas as pd
import seaborn as sns
import numpy as np
import re
import sys
import warnings
import argparse
import matplotlib.ticker as ticker
from collections import defaultdict
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
try:
    from intronCompression import intronCompression
except ModuleNotFoundError:
    from isoespy.intronCompression import intronCompression

def parse_metadata(meta_data, gene):
    ''' metadata file をパースして必要な情報を取得しておく '''
    with open(meta_data) as f:
        lines = f.readlines()
    current_section = None
    sample_meta = dict()
    config_meta = dict()
    gtf_meta = dict()
    colors = dict()
    query = {"gene": [None, None], "tx": [None, None]}
    feature_meta = dict()
    for line in lines:
        line = line.strip()
        if line.startswith("[") and line.endswith("]"):
            current_section = line[1:-1]

        elif current_section == "sample":
            if line.startswith("#"):
                ctrl_trgt = line.lstrip("#")
                ctrl_trgt = [i.strip() for i in ctrl_trgt.split(",")]
                ctrl_trgt_table = dict()
                for pair in ctrl_trgt:
                    CT, var = pair.split("=")
                    ctrl_trgt_table[var.strip()] = CT.strip().lower()
            else:
               if line:
                   if "\t" in line:
                       i = line.split("\t")
                       sample_meta[i[0]] = i[1]
                   else:
                       i = line.split(",")
                       sample_meta[i[0]] = i[1]

        elif current_section == "config":
            if line.startswith("#"):
                key, value = line.split("=")
                key = key.lstrip("#")
                key = key.strip()
                value = value.strip()
                if value == "":
                    continue
                if key == 'qval':
                    value = float(value)
                    config_meta[key] = value
                elif key == 'order':
                    value = [i.strip() for i in value.split(',')]
                    config_meta[key] = value
                elif key == 'colors':
                    value1, value2 = value.split(":")[0].strip(), value.split(":")[1].strip().split(",")
                    value2 = [i.strip() for i in value2]
                    for tx in value2:
                        colors[tx] = value1 

        elif current_section == "gtf":
            if line.startswith("#"):
                key, value = line.split("=")
                key = key.lstrip("#")
                key = key.strip()
                value = value.strip()
                if value == "":
                    continue
                gtf_meta[key] = value

        elif current_section == "query":
            if line.startswith("#"):
                key, value = line.split("=")
                key = key.lstrip("#")
                key = key.strip()
                value = value.strip()
                if value == "":
                    continue
                if key == 'gene':
                    query[key] = [value, gene]   # e.g., ["gene_symbol", "AGMO"]
                elif key == 'tx':
                    value = value.split(",")
                    value = set([i.strip() for i in value])
                    query[key][1] = value   # e.g., [*, ("tx1", "tx2", "tx3")]

        elif current_section == "features":
            if line.startswith("#"):
                line = line.lstrip("#").strip()
                if not ":" in line:
                    feature_class = line
                    feature_id = None
                elif line.split(":")[1].strip() == "":
                    feature_class = line.split(":")[0].strip()
                    feature_id = None
                else:
                    feature_class, feature_id = line.split(":")
                    feature_class, feature_id = feature_class.strip(), feature_id.strip()
                feature_meta[feature_class] = feature_id

    query['tx'][0] = gtf_meta['transcript_id']
    if query['tx'][1] == set():
        query['tx'][1] = None
    return sample_meta, config_meta, ctrl_trgt_table, gtf_meta, colors, query, feature_meta


def get_isoform_model(gtf_file, gtf_meta, colors, query):
    ''' アイソフォームモデルの作成
　　　　exon座標/CDS座標・ストランド・染色体を取得する
        この関数はexonとCDS 専用'''
    # 初期値
    transcripts = {}        # for exons
    transcripts_CDS = {}    # for CDS

    # target_gene と target_tx を定義
    target_gene = query["gene"][1]
    target_tx = query["tx"][1]

    # open GTF
    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            chrom = chrom.replace("chr","")
            if feature == gtf_meta['exon']:
                attr_dict = {match.group(1): match.group(2) for match in re.finditer(r'(\S+)\s+"([^"]+)"', attributes)}   # 属性フィールドのパース
                transcript_id = attr_dict.get(gtf_meta['transcript_id'])

                # 対象遺伝子以外ならパス
                line_gene = attr_dict.get(query["gene"][0])
                line_tx = attr_dict.get(query["tx"][0])
                if line_gene != target_gene:
                    continue
                if target_tx != None:
                    if not line_tx in target_tx:
                        continue

                # エクソンを加える
                if transcript_id not in transcripts:
                    if strand == "+":
                        transcripts[transcript_id] = [[], 1, chrom]
                    else:
                        transcripts[transcript_id] = [[], -1, chrom]
                transcripts[transcript_id][0].append((int(start), int(end)))

            if feature == gtf_meta['cds']:
                attr_dict = {match.group(1): match.group(2) for match in re.finditer(r'(\S+)\s+"([^"]+)"', attributes)}
                transcript_id = attr_dict.get(gtf_meta['transcript_id'])

                # 対象遺伝子以外ならパス
                line_gene = attr_dict.get(query["gene"][0])
                line_tx = attr_dict.get(query["tx"][0])
                if line_gene != target_gene:
                    continue
                if target_tx != None:
                    if not line_tx in target_tx:
                        continue

                # CDSを加える
                if transcript_id not in transcripts_CDS:
                    if strand == "+":
                        transcripts_CDS[transcript_id] = [[], 1, chrom]
                    else:
                        transcripts_CDS[transcript_id] = [[], -1, chrom]
                transcripts_CDS[transcript_id][0].append((int(start), int(end)))

    # 座標順に並び替え
    for tx in transcripts:
        transcripts[tx][0] = sorted(transcripts[tx][0], key=lambda x:x[0])
    for tx in transcripts_CDS:
        transcripts_CDS[tx]

    # colorsを補完
    for tx in transcripts:
        if not tx in colors:
            colors[tx] = "#B3C8CF"
    return transcripts, transcripts_CDS, colors


def get_feature_model(gtf_file, gtf_meta, query, feature_meta):
    ''' feature モデルの取得
        ff_d = {tx1:
                    {"protein_domain":
                                      {id1: [(O,O),(O,O),...(O,O)], id2: [], id3: []}
                     "signal_peptide":
                                      {                                             }
                     ...
                    }
                tx2:
                    {
                    ...
                    }
               }                                                                         '''

    # 初期設定
    ff_d = dict()

    # target_gene と target_tx を定義
    target_gene = query["gene"][1]
    target_tx = query["tx"][1]

    # open GTF
    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            chrom = chrom.replace("chr","")
            if feature in feature_meta:
                attr_dict = {match.group(1): match.group(2) for match in re.finditer(r'(\S+)\s+"([^"]+)"', attributes)}   # 属性フィールドのパース
                transcript_id = attr_dict.get(gtf_meta['transcript_id'])

                # 対象遺伝子以外ならパス
                line_gene = attr_dict.get(query["gene"][0])
                line_tx = attr_dict.get(query["tx"][0])
                if line_gene != target_gene:
                    continue
                if target_tx != None:
                    if not line_tx in target_tx:
                        continue

                # feature 処理
                # tx名が存在しなければ作成
                if transcript_id not in ff_d:
                    ff_d[transcript_id] = {feat: {} for feat in feature_meta}
                # 名無featureの場合
                if feature_meta[feature] == None:
                    if ff_d[transcript_id][feature] == {}:
                        ff_d[transcript_id][feature]["NONAME"] = []
                    ff_d[transcript_id][feature]["NONAME"].append((int(start), int(end)))
                # 名有featureの場合
                else:
                    feat_id = attr_dict.get(feature_meta[feature])
                    if not feat_id in ff_d[transcript_id][feature]:
                        ff_d[transcript_id][feature][feat_id] = []
                    ff_d[transcript_id][feature][feat_id].append((int(start), int(end)))

    # 座標順に並び替え
    for tx in ff_d:
        for feat in feature_meta:
            for ind in ff_d[tx][feat]:
                ff_d[tx][feat][ind] = sorted(ff_d[tx][feat][ind], key=lambda x:x[0])

    return ff_d


def formatting_isoform_model(transcripts_data, transcripts, annot):
    if annot == "exon":
        for transcript_id, exons in transcripts.items():
            start = min([i[0] for i in exons[0]])
            end = max(i[1] for i in exons[0])
            transcripts_data.append({'id': transcript_id, 'exons': exons[0], 'strand': exons[1], 'seq_region_name': exons[2], 'start': start, 'end': end})
    elif annot == "cds":
        for i in range(len(transcripts_data)):
            isomodel = transcripts_data[i]
            tx_id = isomodel['id']
            if tx_id in transcripts:
                # CDSを持つアイソフォーム
                transcripts_data[i]['cds'] = transcripts[tx_id][0]
            else:
                # CDSを持たないアイソフォーム
                transcripts_data[i]['cds'] = []
    return transcripts_data


def get_nonGTF_annotations(annotation_file, gene):
    """	#transcript_id		gene_id	feature_name	feature_type	status
	ENST00000342526.8	AGMO	NMD		binary		0
	AGMO_novel_tx_1		AGMO	NMD		binary		1
	AGMO_novel_tx_2		AGMO	NMD		binary		0
	ENST00000342526.8	AGMO	exp_T		continuous	100
	AGMO_novel_tx_1		AGMO	exp_T		continuous	10
	AGMO_novel_tx_2		AGMO	exp_T		continuous	1
	ENST00000342526.8	AGMO	localization	categorical	A
	AGMO_novel_tx_1		AGMO	localization	categorical	B
	AGMO_novel_tx_2		AGMO	localization	categorical	C
	AGMO_novel_tx_3		AGMO	localization	categorical	C
    """
    if annotation_file == None:
        return {}

    nongtf_d = dict()		# tx_id が key
    with open(annotation_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or line == '':
                continue
            tx_id, gene_id, feat_name, feat_type, status = line.split("\t")
            # 解析遺伝子のみ対象
            if gene != gene_id:
                continue
            # データ加える
            if not tx_id in nongtf_d:
                nongtf_d[tx_id] = {}
            nongtf_d[tx_id][feat_name] = [feat_type, status]

    features_d = dict()		# feature_name が key
    for tx_id in nongtf_d:
        for feat in nongtf_d[tx_id]:
            feat_type = nongtf_d[tx_id][feat][0]
            if not feat in features_d:
                features_d[feat] = [feat_type, {}]
    for feat in features_d:
        features_d[feat][1] = {tx_id: None for tx_id in nongtf_d}   #初期値
    for tx_id in nongtf_d:
        for feat in nongtf_d[tx_id]:
            status = nongtf_d[tx_id][feat][1]
            features_d[feat][1][tx_id] = status		#statusが存在する場合はNoneを値で上書き

    return features_d


def prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff):
    if ci == None:
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    else:
        ax1.set_xticks([x_min_eff, x_max_eff])
        ax1.set_xticklabels([str(x_min), str(x_max)])
    return ax1


def coordinate_adjustment(feat_names, Y_LOW, Y_UPP):
    # featureごとのy範囲を返す
    """ coord_d = {"protein_coding": (y_low, y_upp),
                   "signal_peptide": (y_low, y_upp),
                   ... }                               """
    # out var
    out_d = dict()
    if feat_names == []:
        return out_d
    # divide
    width = (Y_UPP - Y_LOW) / len(feat_names)
    for i in range(len(feat_names)):
        out_d[feat_names[i]] = (Y_UPP-width*(i+1), Y_UPP-width*i)
    return out_d


def darken_color(color, factor=0.4):
    """
    RGBカラーを少し暗くする関数。
    factor: 0.0～1.0の範囲で調整（小さいほど暗くなる）
    """
    return tuple(max(0, c * factor) for c in color)


def overlaped_feature_indivisuals(DICT):
    # 重なり判定関数
    def overlaps(r1, r2):
        return r1[0] <= r2[1] and r2[0] <= r1[1]

    # 1. グラフ構築：重なりがあるペアでエッジを作る
    graph = defaultdict(set)
    keys = list(DICT.keys())

    for i in range(len(keys)):
        for j in range(i+1, len(keys)):
            k1, k2 = keys[i], keys[j]
            for r1 in DICT[k1]:
                for r2 in DICT[k2]:
                    if overlaps(r1, r2):
                        graph[k1].add(k2)
                        graph[k2].add(k1)
                        break
                else:
                    continue
                break

    # 2. DFSで連結成分（グループ）を作る
    visited = set()
    groups = []

    def dfs(node, group):
        visited.add(node)
        group.append(node)
        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor, group)

    for k in keys:
        if k not in visited:
            group = []
            dfs(k, group)
            groups.append(group)

    return groups


def prepare_ff_ax1(transcripts_data, model_features, ax1, gene_name, colors):
    # AX1: functional features を視覚化
    # 準備
    y_positions = []							   # 転写産物名の記入のため
    MIN = min(transcripts_data, key=lambda x: x['start'])['start']	   # ターゲット遺伝子全転写産物の最左端
    MAX = max(transcripts_data, key=lambda x: x['end'])['end']		   # ターゲット遺伝子全転写産物の最右端
    transcripts_data = transcripts_data[::-1]				   # 作図の順番の都合上
    if model_features == {}:
        feat_names = []
        feat_colors = []
    else:
        feat_names = list(model_features[next(iter(model_features))].keys())   # metadata 指定のfeature名のリスト
        feat_colors = sns.color_palette("husl", len(feat_names))

    for i, transcript_data in enumerate(transcripts_data):
        y_positions.append(i)
        start = transcript_data['start']
        end = transcript_data['end']
        tx_name = transcript_data['id']

        # 1. 背骨 を描画 (エクソンの背面に灰色の矢印付きの線分を描画)
        strand = transcript_data['strand']
        arrow_direction = "right" if strand == 1 else "left"
        y_pos = i
        interval = (MAX-MIN)//50
        x_positions = np.arange(start, end, interval)
        x_positions = x_positions[1:-1]

        ax1.annotate('', xy=(end, y_pos), xytext=(start, y_pos), arrowprops=dict(arrowstyle="-", color='gray', lw=1))
        if arrow_direction == "right":
            for x in x_positions:
                ax1.scatter(x, y_pos, marker=">", color="gray", s=10)
        else:
            for x in x_positions:
                ax1.scatter(x, y_pos, marker="<", color="gray", s=10)

        # 2. exon を描画（少し濃い水色で描画）
        for exon in transcript_data['exons']:
            exon_start = exon[0]
            exon_end = exon[1]
            # CDS以外の部分を新しい色で描画
            ax1.add_patch(patches.Rectangle((exon_start, i - 0.04), exon_end - exon_start, 0.08, color=colors[tx_name]))

        # 3. CDS を描画
        for cds in transcript_data['cds']:
            cds_start = cds[0]
            cds_end = cds[1]
            ax1.add_patch(patches.Rectangle((cds_start, i - 0.08), cds_end - cds_start, 0.16, color=colors[tx_name]))

        # 4. feature を描画 (ボックス)
        if not tx_name in model_features:
            # featureが1つも存在しないケース
            continue
        else:
            # featureが少なくとも1つ存在するケース
            Y_LOW = i - 0.8
            Y_UPP = i - 0.11
            coord_d = coordinate_adjustment(feat_names, Y_LOW, Y_UPP)
            for i, feat in enumerate(feat_names):
                # for features of an isoform
                y_low, y_upp = coord_d[feat][0], coord_d[feat][1]
                # 範囲の降順: featureの重なりによる見えにくさを考慮するため大きい順にプロットする
                plot_order = sorted(model_features[tx_name][feat], key=lambda k: sum(e - s for s, e in model_features[tx_name][feat][k]), reverse=True)
                # 重なりのあるfeatureをグルーピング
                plot_groups = overlaped_feature_indivisuals(model_features[tx_name][feat])
                plot_groups_cnt = [-1 for i in range(len(plot_groups))]
                # 実際のプロット
                for ind in plot_order:
                    # for individuals of a feature
                    # プロットするy座標
                    for g in range(len(plot_groups)):
                        if ind in plot_groups[g]:
                            plot_groups_cnt[g] += 1
                            idx = g   # 重なる場合に少し下にずらすためのカウント
                            break
                    y_low2 = y_low - 0.009*plot_groups_cnt[idx]
                    y_upp2 = y_upp - 0.009*plot_groups_cnt[idx]
                    crush_rate = 0.83
                    # step 1: スプライシングパターンに従って色塗り
                    for k in model_features[tx_name][feat][ind]:
                        start_k, end_k = k[0], k[1]
                        ax1.add_patch(patches.Rectangle((start_k, y_upp2-(y_upp2-y_low2)*crush_rate), end_k - start_k, (y_upp2 - y_low2)*crush_rate, color=feat_colors[i]))
                    # step 2: 要素ごとに区画化
                    leftmost_x, rightmost_x = model_features[tx_name][feat][ind][0][0], model_features[tx_name][feat][ind][-1][1]
                    ax1.add_patch(patches.Rectangle((leftmost_x, y_upp2-(y_upp2-y_low2)*crush_rate), rightmost_x - leftmost_x, (y_upp2 - y_low2)*crush_rate, linewidth=0.5, edgecolor=darken_color(feat_colors[i]), facecolor='none'))


    # 軸と図の設定
    space = int((MAX-MIN)/20)
    ax1.set_xlim(MIN - 2*space, MAX + space)
    ax1.set_ylim(-1.0, len(transcripts_data) - 0.5)
    ax1.set_xlabel(f'Chr{transcript_data["seq_region_name"]}')
    ax1.set_title(gene_name)
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    ax1.set_yticks(y_positions)
    transcripts = [i['id'] for i in transcripts_data]
    ax1.set_yticklabels(transcripts)

    # feature legend
    N = len(transcripts_data) - 1
    coord = coordinate_adjustment(feat_names, Y_LOW=N-0.8, Y_UPP=N-0.11)
    for i, feat in enumerate(feat_names):
        y_low, y_upp = coord[feat][0], coord[feat][1]
        ax1.text(MIN-1.8*space, (y_upp+y_low)/2, feat, fontsize=9, color=darken_color(feat_colors[i], factor=0.6), ha='left', va='center')
    return ax1


def prepare_ff_axk(ax_l, transcripts_data, annotation_dataB):
    def Color_binary(A):
        # Binaryの色マッピング
        binary_colors = {"0": (0.875, 0.9, 0.875, 0.85), "1": (0.4, 0.7, 0.4, 0.9), None: 'white'}
        # Aの値に対応する色を取得
        A_colors = {key: binary_colors[value] for key, value in A.items()}
        return A_colors

    def Color_categorical(A):
        # Categoricalの色マッピング
        unique_categories = sorted(set(v for v in A.values() if v is not None))
        cmap = plt.get_cmap("Pastel1", len(unique_categories))
        category_colors = {cat: cmap(i) for i, cat in enumerate(unique_categories)}
        category_colors[None] = 'white'
        A_colors = {key: category_colors[value] for key, value in A.items()}
        return A_colors

    def is_float(v):
        try:
            float(v)
            return True
        except (ValueError, TypeError):
            return False

    def Color_continuous(A):
        # Continuousの色マッピング
        valid_values = [float(v) for v in A.values() if v is not None and is_float(v)]		#[v for v in A.values() if v is not None]
        vmin, vmax = min(valid_values), max(valid_values)
        cmap = plt.get_cmap("Reds")
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        A_colors = {key: cmap(norm(float(value))) if value is not None and is_float(value) else 'white' for key, value in A.items()}
        return A_colors

    # 準備(データ補完)
    for transcript_data in transcripts_data:
        tx_name = transcript_data['id']
        for feat in annotation_dataB:
            if not tx_name in annotation_dataB[feat][1]:
                annotation_dataB[feat][1][tx_name] = None

    # axに追加
    N = len(transcripts_data)
    features_l = list(annotation_dataB.keys())
    for i, transcript_data in enumerate(transcripts_data):
        # 横向きに処理(isoformごとに処理)
        tx_name = transcript_data['id']
        # 作図
        for j, feat in enumerate(features_l):
            CATEGORY = annotation_dataB[feat][0]
            VALUES = annotation_dataB[feat][1]
            if CATEGORY == "binary":
                # binary featureの場合のカラーマップ
                colormap = Color_binary(VALUES)
                if VALUES[tx_name] == None:
                    rect = patches.Rectangle((0.3, N-i-0.9-1), 0.4, 1, facecolor=colormap[tx_name], edgecolor='black', hatch='//', linewidth=2)
                else:
                    rect = patches.Rectangle((0.3, N-i-0.9-1), 0.4, 1, facecolor=colormap[tx_name], edgecolor='black', linewidth=2)
                ax_l[j].add_patch(rect)
                # 文字を記入
                x = 0.3 + 0.4 / 2        # 長方形のx中心
                y = N - i - 0.9 - 1 + 1/2 # 長方形のy中心
                ax_l[j].text(x, y, VALUES[tx_name], ha='center', va='center', fontsize=10)
            elif CATEGORY == "categorical":
                # categorical featureの場合のカラーマップ
                colormap = Color_categorical(VALUES)
                if VALUES[tx_name] == None:
                    rect = patches.Rectangle((0.3, N-i-0.9-1), 0.4, 1, facecolor=colormap[tx_name], edgecolor='black', hatch='//', linewidth=2)
                else:
                    rect = patches.Rectangle((0.3, N-i-0.9-1), 0.4, 1, facecolor=colormap[tx_name], edgecolor='black', linewidth=2)
                ax_l[j].add_patch(rect)
                # 文字を記入
                x = 0.3 + 0.4 / 2        # 長方形のx中心
                y = N - i - 0.9 - 1 + 1/2 # 長方形のy中心
                ax_l[j].text(x, y, VALUES[tx_name], ha='center', va='center', rotation='vertical', fontsize=10)
            elif CATEGORY == "continuous":
                # continuous featureの場合のカラーマップ
                colormap = Color_continuous(VALUES)
                if VALUES[tx_name] == None:
                    rect = patches.Rectangle((0.3, N-i-0.9-1), 0.4, 1, facecolor=colormap[tx_name], edgecolor='black', hatch='//', linewidth=2)
                else:
                    rect = patches.Rectangle((0.3, N-i-0.9-1), 0.4, 1, facecolor=colormap[tx_name], edgecolor='black', linewidth=2)
                ax_l[j].add_patch(rect)
                # 文字を記入
                x = 0.3 + 0.4 / 2        # 長方形のx中心
                y = N - i - 0.9 - 1 + 1/2 # 長方形のy中心
                ax_l[j].text(x, y, VALUES[tx_name], ha='center', va='center', rotation='vertical', fontsize=10)

    for j, feat in enumerate(features_l):
        ax_l[j].set_title(feat)
        #ax_l[j].text(0.5, N-1, feat, transform=ax_l[j].transAxes, ha='center', va='bottom', fontsize=9)

    # 位置やtickの設定
    for ax in ax_l:
        ax.set_xlim(0, 1)
        ax.set_ylim(-1.0, N - 0.5)
        ax.axis('off')
    return ax_l


def reorder(transcripts_data, meta_data):
    if not 'order' in meta_data:
        return transcripts_data
    else:
        tmp = []
        order = meta_data['order']
        for tx in order:
            for i in transcripts_data:
                if i['id'] == tx:
                    tmp.append(i)
                    break
        return tmp


def plot_isoespy_ff(transcripts_data, model_features, config_meta, gene_name, meta_data, colors, ci, x_min, x_max, x_min_eff, x_max_eff, annotation_dataB):
    # 並び替え
    transcripts_data = reorder(transcripts_data, meta_data)

    # プロットする場所を作成
    N = len(annotation_dataB.keys())    # non-positional featuresの種類
    if N == 0:
        fig, ax1 = plt.subplots(ncols=1, figsize=(20, 8), gridspec_kw={'width_ratios': [20] + [1] * N})
    else:
        fig, axes = plt.subplots(ncols=N+1, figsize=(20, 8), gridspec_kw={'width_ratios': [20] + [1] * N})
        ax1, *ax_l = axes

    # prep for ax1
    ax1 = prepare_ff_ax1(transcripts_data, model_features, ax1, gene_name, colors)
    ax1 = prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff)

    # prep for axk
    if N != 0:
        ax_l = prepare_ff_axk(ax_l, transcripts_data, annotation_dataB)

    plt.subplots_adjust(wspace=0.02)
    plt.show()


def moved_data_for_exons_cds(main_data, ci):
    """ intron compression 処理: exons, cds, start, end を補正 """
    # preparation
    model = {}
    for tx_data in main_data:
        tx_name = tx_data['id']
        model[tx_name+"_exons"] = tx_data['exons']
        model[tx_name+"_cds"] = tx_data['cds']

    # intron compression
    model_compressed = intronCompression(model, ci)

    # start, end
    startend_d = {}
    for tx_data in main_data:
        tx_name = tx_data['id']
        start = model_compressed[tx_name+"_exons"][0][0]
        end = model_compressed[tx_name+"_exons"][-1][1]
        startend_d[tx_name] = {"start": start, "end": end}

    # main_dataの更新
    for tx_data in main_data:
        tx_name = tx_data['id']
        tx_data['exons'] = model_compressed[tx_name+"_exons"]
        tx_data['cds'] = model_compressed[tx_name+"_cds"]
        tx_data['start'] = startend_d[tx_name]["start"]
        tx_data['end'] = startend_d[tx_name]["end"]

    return main_data


def moved_data_for_features(model_features, main_data, ci):
    """ intron compression 処理: features の座標を補正 """
    # preparation
    model = {}
    for tx_data in main_data:
        tx_name = tx_data['id']
        model[tx_name+"_exons"] = tx_data['exons']
        model[tx_name+"_cds"] = tx_data['cds']
    for tx_name in model_features:
        for feat_name in model_features[tx_name]:
            for id_name in model_features[tx_name][feat_name]:
                model[tx_name+"_"+feat_name+"_"+id_name] = model_features[tx_name][feat_name][id_name]

    # intron compression
    model_compressed = intronCompression(model, ci)

    # model_featuresの更新
    for tx_name in model_features:
        for feat_name in model_features[tx_name]:
            for id_name in model_features[tx_name][feat_name]:
                model_features[tx_name][feat_name][id_name] = model_compressed[tx_name+"_"+feat_name+"_"+id_name]

    return model_features


def isoespy_ff(gene, gtf_data, meta_data, ci, annotation_file):
    # metadata
    sample_meta, config_meta, ctrl_trgt_table, gtf_meta, colors, query, feature_meta = parse_metadata(meta_data, gene)

    # isoform model
    model_exon, model_cds, colors = get_isoform_model(gtf_data, gtf_meta, colors, query)

    # feature model
    model_features = get_feature_model(gtf_data, gtf_meta, query, feature_meta)

    # formatted isoform model
    transcripts_main_data = list()
    transcripts_main_data = formatting_isoform_model(transcripts_main_data, model_exon, annot="exon")
    transcripts_main_data = formatting_isoform_model(transcripts_main_data, model_cds, annot="cds")
    #print("transcripts_main_data", transcripts_main_data)
    x_min = min(transcripts_main_data, key=lambda x: x['start'])['start']
    x_max = max(transcripts_main_data, key=lambda x: x['end'])['end']

    if ci != None:
        # イントロン圧縮する場合
        copied_data = copy.deepcopy(transcripts_main_data)
        transcripts_main_data = moved_data_for_exons_cds(transcripts_main_data, ci)
        # ff()では加えてfeaturesの補正も必要
        model_features = moved_data_for_features(model_features, copied_data, ci)
    x_min_eff = min(transcripts_main_data, key=lambda x: x['start'])['start']
    x_max_eff = max(transcripts_main_data, key=lambda x: x['end'])['end']

    # non-GTF annotation features
    annotation_dataB = get_nonGTF_annotations(annotation_file, gene)

    # plot
    plot_isoespy_ff(transcripts_main_data, model_features, config_meta, gene, config_meta, colors, ci, x_min, x_max, x_min_eff, x_max_eff, annotation_dataB)


def process_ci(ci):
    if ci is None:
        return None  # デフォルトのまま

    if isinstance(ci, str):
        try:
            ci = float(ci)
        except ValueError:
            warnings.warn(f"Warning: ci should be a float/int, but received string '{ci}' that cannot be converted.")
            return None  # 変換できない場合は None にする
    else:
        try:
            ci = float(ci)  # float に変換できるなら変換
        except (ValueError, TypeError):
            warnings.warn(f"Warning: ci should be a float/int, but received '{ci}' that cannot be converted.")
            return None  # 変換できない場合は None にする

    return ci  # 変換された float 値を返す


def main(args=None):
    # parser
    parser = argparse.ArgumentParser(description='isoespy_ff()')
    parser.add_argument('-gene', '--gene_name', required=True, type=str, default=None, help='Gene name')
    parser.add_argument('-gtf', '--gtf_data', required=True, type=str, default=None, help='GTF file')
    parser.add_argument('-meta', '--meta_data', required=True, type=str, default=None, help='metadata')
    parser.add_argument('-ci', '--compress_introns', default=None, help='intron compression parameter')
    parser.add_argument('-a', '--annotation', type=str, default=None, help='non-gtf transcript annotation file')

    args = parser.parse_args()
    gene = args.gene_name
    gtf_data = args.gtf_data
    meta_data = args.meta_data
    ci = args.compress_introns    # None(デフォルト)ならcompressionしない。float
    ci = process_ci(ci)
    annotation_data = args.annotation

    # isoespy_ff()
    isoespy_ff(gene, gtf_data, meta_data, ci, annotation_data)


if __name__ == '__main__':
    main()
