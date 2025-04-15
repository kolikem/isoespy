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
                    query["gene"] = [value, gene]   # e.g., ["gene_symbol", "AGMO"]
                elif key == 'tx':
                    value = value.split(",")
                    value = set([i.strip() for i in value])
                    query["tx"][1] = value   # e.g., [*, ("tx1", "tx2", "tx3")]

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


def get_expression_data(expression_file, transcripts_data, ctrl_trgt_table, sample_meta):
    df = pd.read_csv(expression_file, sep="\t", header=0)
    targets = [i['id'] for i in transcripts_data]
    filtered_df = df[df.iloc[:, 0].isin(targets)]
    GROUP_0, GROUP_1 = dict(), dict()
    for index, row in filtered_df.iterrows():
        tx_id = row.iloc[0]
        group_0, group_1 = [], []
        for sample in filtered_df.columns[1:]:
            expression = row[sample]
            group_symbol = sample_meta[sample]
            if ctrl_trgt_table[group_symbol] == 'control':
                group_0.append(expression)
            elif ctrl_trgt_table[group_symbol] == 'target':
                group_1.append(expression)

        GROUP_0[tx_id] = group_0
        GROUP_1[tx_id] = group_1

    # integrate expression data to main data
    for i in range(len(transcripts_data)):
        transcripts_data[i]['group0_exp'] = GROUP_0[transcripts_data[i]['id']]
        transcripts_data[i]['group1_exp'] = GROUP_1[transcripts_data[i]['id']]

    return transcripts_data


def get_det_data(det_file, transcripts_data):
    det_data = dict()
    with open(det_file) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                if "\t" in line:
                    col = line.split("\t")
                else:
                    col = line.split(",")
                tx_id = col[0]
                logFC = float(col[1])
                qval = float(col[2])
                det_data[tx_id] = [logFC, qval]

    # integrate det data to main data
    for i in range(len(transcripts_data)):
        transcripts_data[i]['det'] = det_data[transcripts_data[i]['id']]

    return transcripts_data


def prepare_de_ax1(transcripts_data, ax1, gene_name, colors):
    # AX1: 各アイソフォームを視覚化
    y_positions = []
    MIN = min(transcripts_data, key=lambda x: x['start'])['start']
    MAX = max(transcripts_data, key=lambda x: x['end'])['end']
    transcripts_data = transcripts_data[::-1]
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
            ax1.add_patch(patches.Rectangle((exon_start, i - 0.1), exon_end - exon_start, 0.2, color=colors[tx_name]))

        # 3. CDS を描画
        for cds in transcript_data['cds']:
            cds_start = cds[0]
            cds_end = cds[1]
            ax1.add_patch(patches.Rectangle((cds_start, i - 0.2), cds_end - cds_start, 0.4, color=colors[tx_name]))

    # 軸と図の設定
    space = int((MAX-MIN)/20)
    ax1.set_xlim(MIN - space, MAX + space)
    ax1.set_ylim(-0.5, len(transcripts_data) - 0.5)
    ax1.set_xlabel(f'Chr{transcript_data["seq_region_name"]}')
    ax1.set_title(gene_name)
    #ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    ax1.set_yticks(y_positions)
    transcripts = [i['id'] for i in transcripts_data]
    ax1.set_yticklabels(transcripts)
    return ax1


def prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff):
    if ci == None:
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    else:
        ax1.set_xticks([x_min_eff, x_max_eff])
        ax1.set_xticklabels([str(x_min), str(x_max)])
    return ax1


def prepare_de_ax2(transcripts_data, ax2, config_meta, outliers):
    # AX2: 箱ひげ図
    transcripts_data = transcripts_data[::-1]
    transcripts = [i['id'] for i in transcripts_data]
    labels = []
    data = []
    for i in range(len(transcripts_data)):
        tx_id = transcripts_data[i]['id']
        data.append(transcripts_data[i]['group1_exp'])
        data.append(transcripts_data[i]['group0_exp'])
        labels.append(f'{tx_id}_1')
        labels.append(f'{tx_id}_0')

    # プロットの位置調整
    positions = []
    for i, name in enumerate(transcripts):
    #    positions.extend([i*3+1, i*3+1.6])
         positions.extend([i-0.2, i+0.2])

    # 箱ひげ図をプロット
    boxplot = ax2.boxplot(data, positions=positions, vert=False, patch_artist=True, widths=0.2, showfliers=outliers)

    # 色を設定 (奇数: 薄い, 偶数: 濃い)
    colors = ['lightgray' if idx % 2 == 1 else 'gray' for idx in range(len(data))]
    for patch, color in zip(boxplot['boxes'], colors):
        patch.set_facecolor(color)

    # 詳細設定
    for median in boxplot['medians']:
        median.set_color('black')  # 色を黒に設定

    ax2.set_ylim(-0.5, len(transcripts_data) - 0.5)
    ax2_xlabel = "Expression level"
    #if 'ax2_xlabel' in config_meta:
    #    ax2_xlabel = config_meta['ax2_xlabel']
    ax2.set_yticklabels([])
    ax2.tick_params(axis='y', length=0)
    ax2.set_xlabel(ax2_xlabel)
    return ax2


def prepare_de_ax3(transcripts_data, config_meta, ax3):
    # AX3: 色塗り
    def color_mapping(nlogqval):
        if nlogqval > 20:
            return 1
        else:
            return 0.3 + nlogqval*0.035
    DET = list()
    for i in range(len(transcripts_data)):
        DET.append([transcripts_data[i]['id']] + transcripts_data[i]['det'])
    qval_threshold = config_meta['qval']
    # カラーマップを設定
    warm_cmap = plt.cm.Reds    # 暖色系
    cold_cmap = plt.cm.Blues   # 寒色系
    # 最大値のスケーリングを計算（色の濃さを調整するため）
    max_value = max(-np.log10(v[2]) for v in DET)
    for i, val in enumerate(DET):
        logFC = DET[i][1]
        qval = DET[i][2]
        # 色の設定
        if qval < qval_threshold:
            if logFC > 0:
                color = warm_cmap(color_mapping(-np.log10(qval)))
            else:
                color = cold_cmap(color_mapping(-np.log10(qval)))
        else:
                color = (0.5, 0.5, 0.5, 0.5)  # Gray
        # 長方形の作図
        rect = patches.Rectangle((0, len(DET) - i - 1), 0.1, 1, facecolor=color, edgecolor='black', linewidth=2)
        ax3.add_patch(rect)

    # 凡例の作成
    reds_colors = [plt.cm.Reds(x) for x in np.linspace(0.3, 1, 128)]
    blues_colors = [plt.cm.Blues(x) for x in np.linspace(0.3, 1, 128)]
    combined_colors = blues_colors[::-1] + reds_colors
    combined_cmap = LinearSegmentedColormap.from_list("CombinedMap", combined_colors)
    norm = Normalize(vmin=0, vmax=1)
    sm = ScalarMappable(cmap=combined_cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax3, orientation='vertical', fraction=0.05, pad=0.01, shrink=0.8)
    cbar.set_label('-log10(q)', rotation=90)
    cbar.set_ticks(np.linspace(0, 1, 5))
    tick_labels = [">20", "10", "0", "10", ">20"]
    cbar.ax.set_yticklabels(tick_labels)

    # Set limits and remove axes for a cleaner look
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, len(DET))
    ax3.axis('off')
    return ax3


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


def plot_isoespy_de(transcripts_data, config_meta, gene_name, meta_data, colors, ci, x_min, x_max, x_min_eff, x_max_eff, outliers):
    # 並び替え
    transcripts_data = reorder(transcripts_data, meta_data)

    # プロットする場所を作成
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(20, 5), gridspec_kw={'width_ratios': [10, 6, 2]})

    # prep for ax1
    ax1 = prepare_de_ax1(transcripts_data, ax1, gene_name, colors)
    ax1 = prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff)

    # prep for ax2
    ax2 = prepare_de_ax2(transcripts_data, ax2, config_meta, outliers)

    # prep for ax3
    ax3 = prepare_de_ax3(transcripts_data, config_meta, ax3)

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


def isoespy_de(gene, gtf_data, expression_data, det_data, meta_data, ci, outliers):
    # metadata
    sample_meta, config_meta, ctrl_trgt_table, gtf_meta, colors, query, feature_meta = parse_metadata(meta_data, gene)

    # raw isoform model
    model_exon, model_cds, colors = get_isoform_model(gtf_data, gtf_meta, colors, query)

    # formatted isoform model
    transcripts_main_data = list()
    transcripts_main_data = formatting_isoform_model(transcripts_main_data, model_exon, annot="exon")
    transcripts_main_data = formatting_isoform_model(transcripts_main_data, model_cds, annot="cds")
    x_min = min(transcripts_main_data, key=lambda x: x['start'])['start']
    x_max = max(transcripts_main_data, key=lambda x: x['end'])['end']

    # transcripts_main_data の補正
    if ci != None:
        # イントロン圧縮あり
        transcripts_main_data = moved_data_for_exons_cds(transcripts_main_data, ci)
    x_min_eff = min(transcripts_main_data, key=lambda x: x['start'])['start']
    x_max_eff = max(transcripts_main_data, key=lambda x: x['end'])['end']

    # epxression data
    transcripts_main_data = get_expression_data(expression_data, transcripts_main_data, ctrl_trgt_table, sample_meta)

    # det data
    transcripts_main_data = get_det_data(det_data, transcripts_main_data)

    # plot
    plot_isoespy_de(transcripts_main_data, config_meta, gene, config_meta, colors, ci, x_min, x_max, x_min_eff, x_max_eff, outliers)


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
    parser = argparse.ArgumentParser(description='isoespy_de()')
    parser.add_argument('-gene', '--gene_name', required=True, type=str, default=None, help='Gene name')
    parser.add_argument('-gtf', '--gtf_data', required=True, type=str, default=None, help='GTF file')
    parser.add_argument('-exp', '--expression_data', required=True, type=str, default=None, help='Expression data')
    parser.add_argument('-det', '--det_data', required=True, type=str, default=None, help='Differential expression data')
    parser.add_argument('-meta', '--meta_data', required=True, type=str, default=None, help='metadata')
    parser.add_argument('-ci', '--compress_introns', default=None, help='intron compression parameter')
    parser.add_argument('--show_outliers', action="store_true", help="Show outliers in the boxplot (default: hide)")

    args = parser.parse_args()
    gene = args.gene_name
    gtf_data = args.gtf_data
    expression_data = args.expression_data
    det_data = args.det_data
    meta_data = args.meta_data
    ci = args.compress_introns    # None(デフォルト)ならcompressionしない。float
    ci = process_ci(ci)
    outliers = args.show_outliers

    # isoespy_de()
    isoespy_de(gene, gtf_data, expression_data, det_data, meta_data, ci, outliers)


if __name__ == '__main__':
    main()
