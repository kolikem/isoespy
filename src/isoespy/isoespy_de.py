import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as mtransforms
import copy
import pandas as pd
import seaborn as sns  # 追加
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
    try:
        from isoespy.intronCompression import intronCompression
    except ModuleNotFoundError:
        def intronCompression(model, ci):
            return model

# -----------------------------
# metadata parser
# -----------------------------
def parse_metadata(meta_data, gene):
    """Parse metadata file and return required dicts."""
    with open(meta_data) as f:
        lines = f.readlines()

    current_section = None
    sample_meta = dict()
    config_meta = dict()
    gtf_meta = dict()
    query = {"gene": [None, None], "tx": [None, None]}
    feature_meta = dict()

    # Hierarchical colors meta (Matched with ff structure)
    colors_meta = {
        "global": {
            "palette": None,
            "default_line": "gray",
            "default_text": "black",
            "default_tx": "#B3C8CF", # Default Transcript Color
        },
        "transcripts": {},
        "exon": {"color": None},
        "cds": {"color": None},
        # DE specific settings
        "de": {
            "box_group0": "#4daf4a", # Default control color
            "box_group1": "#e41a1c"  # Default target color
        },
    }

    ctrl_trgt_table = dict()

    for line in lines:
        line = line.strip()
        # コメント行の処理 (!で始まる行は無視)
        if not line or line.startswith("!"):
            continue

        if line.startswith("[") and line.endswith("]"):
            current_section = line[1:-1]
            continue

        # パラメータ行の解析 (#key=value または key=value)
        content = line.lstrip("#").strip() if line.startswith("#") else line

        # ---------------- sample ----------------
        if current_section == "sample":
            # sampleセクションは #始まりが target definition, それ以外が sample table
            if line.startswith("#"):
                ctrl_trgt = content.split(",")
                for pair in ctrl_trgt:
                    if "=" in pair:
                        CT, var = pair.split("=")
                        ctrl_trgt_table[var.strip()] = CT.strip().lower()
            else:
                if "\t" in line:
                    i = line.split("\t")
                    sample_meta[i[0]] = i[1]
                elif "," in line:
                    i = line.split(",")
                    sample_meta[i[0]] = i[1]

        # ---------------- config ----------------
        elif current_section == "config":
            if "=" in content:
                key, value = content.split("=", 1)
                key = key.strip()
                value = value.strip()
                if key == "qval":
                    config_meta[key] = float(value)
                elif key == "order":
                    config_meta[key] = [i.strip() for i in value.split(",")]
                elif key == "colors":
                    # legacy style support
                    if ":" in value:
                        val1, val2 = value.split(":", 1)
                        val1 = val1.strip()
                        val2 = [i.strip() for i in val2.split(",")]
                        for tx in val2:
                            if tx not in colors_meta["transcripts"]:
                                colors_meta["transcripts"][tx] = val1

        # ---------------- gtf ----------------
        elif current_section == "gtf":
            if "=" in content:
                key, value = content.split("=", 1)
                gtf_meta[key.strip()] = value.strip()

        # ---------------- query ----------------
        elif current_section == "query":
            if "=" in content:
                key, value = content.split("=", 1)
                key = key.strip()
                value = value.strip()
                if key == "gene":
                    query[key] = [value, gene]
                elif key == "tx":
                    val_set = set([i.strip() for i in value.split(",")])
                    query[key][1] = val_set

        # ---------------- features ----------------
        elif current_section == "features":
            if ":" in content:
                feature_class, feature_id = content.split(":", 1)
                feature_meta[feature_class.strip()] = feature_id.strip() if feature_id.strip() else None
            else:
                feature_meta[content] = None

        # ---------------- colors (General) ----------------
        elif current_section == "colors":
            if "=" in content:
                key, value = content.split("=", 1)
                if key.strip() in colors_meta["global"]:
                    colors_meta["global"][key.strip()] = value.strip()

        # ---------------- colors.transcripts ----------------
        elif current_section == "colors.transcripts":
            if "=" in content:
                key, value = content.split("=", 1)
                colors_meta["transcripts"][key.strip()] = value.strip()

        # ---------------- colors.exon / cds ----------------
        elif current_section == "colors.exon":
            if "=" in content:
                key, value = content.split("=", 1)
                if key.strip() == "color": colors_meta["exon"]["color"] = value.strip()
        
        elif current_section == "colors.cds":
            if "=" in content:
                key, value = content.split("=", 1)
                if key.strip() == "color": colors_meta["cds"]["color"] = value.strip()

        # ---------------- colors.de (Specific) ----------------
        elif current_section == "colors.de":
            if "=" in content:
                key, value = content.split("=", 1)
                key = key.strip()
                if key in colors_meta["de"]:
                    colors_meta["de"][key] = value.strip()

    # finalize query
    query["tx"][0] = gtf_meta.get("transcript_id")
    if query["tx"][1] == set():
        query["tx"][1] = None

    return sample_meta, config_meta, ctrl_trgt_table, gtf_meta, colors_meta, query, feature_meta


def get_isoform_model(gtf_file, gtf_meta, tx_colors, query, colors_meta=None):
    ''' アイソフォームモデルの作成と色の補完 '''
    transcripts = {}        # for exons
    transcripts_CDS = {}    # for CDS
    target_gene = query["gene"][1]
    target_tx = query["tx"][1]

    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            chrom = chrom.replace("chr","")
            
            attr_dict = {match.group(1): match.group(2) for match in re.finditer(r'(\S+)\s+"([^"]+)"', attributes)}
            transcript_id = attr_dict.get(gtf_meta['transcript_id'])
            
            # フィルタリング
            line_gene = attr_dict.get(query["gene"][0])
            line_tx = attr_dict.get(query["tx"][0])
            
            if line_gene != target_gene:
                continue
            if target_tx is not None and line_tx not in target_tx:
                continue

            if feature == gtf_meta['exon']:
                if transcript_id not in transcripts:
                    direction = 1 if strand == "+" else -1
                    transcripts[transcript_id] = [[], direction, chrom]
                transcripts[transcript_id][0].append((int(start), int(end)))

            if feature == gtf_meta.get('cds', 'CDS'):
                if transcript_id not in transcripts_CDS:
                    direction = 1 if strand == "+" else -1
                    transcripts_CDS[transcript_id] = [[], direction, chrom]
                transcripts_CDS[transcript_id][0].append((int(start), int(end)))

    # ソート
    for tx in transcripts:
        transcripts[tx][0] = sorted(transcripts[tx][0], key=lambda x:x[0])
    for tx in transcripts_CDS:
        transcripts_CDS[tx][0] = sorted(transcripts_CDS[tx][0], key=lambda x:x[0])

    # === 色の自動補完 (missing colors) ===
    default_tx = "#B3C8CF"
    palette_name = None
    if colors_meta:
        default_tx = colors_meta["global"].get("default_tx", default_tx)
        palette_name = colors_meta["global"].get("palette", None)

    missing = [tx for tx in transcripts if tx not in tx_colors]
    if missing:
        if palette_name:
            pal = sns.color_palette(palette_name, len(missing))
            for tx, c in zip(missing, pal):
                tx_colors[tx] = mcolors.to_hex(c)
        else:
            # デフォルトパレット (husl)
            pal = sns.color_palette("husl", len(missing))
            for tx, c in zip(missing, pal):
                tx_colors[tx] = mcolors.to_hex(c)
    
    return transcripts, transcripts_CDS, tx_colors


def formatting_isoform_model(transcripts_data, transcripts, annot):
    if annot == "exon":
        for transcript_id, exons in transcripts.items():
            if not exons[0]: continue
            start = min([i[0] for i in exons[0]])
            end = max(i[1] for i in exons[0])
            transcripts_data.append({'id': transcript_id, 'exons': exons[0], 'strand': exons[1], 'seq_region_name': exons[2], 'start': start, 'end': end})
    elif annot == "cds":
        for i in range(len(transcripts_data)):
            isomodel = transcripts_data[i]
            tx_id = isomodel['id']
            if tx_id in transcripts:
                transcripts_data[i]['cds'] = transcripts[tx_id][0]
            else:
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
        # カラム名がサンプル名と仮定
        for sample in filtered_df.columns[1:]:
            if sample not in sample_meta: continue
            expression = row[sample]
            group_symbol = sample_meta[sample]
            
            if group_symbol not in ctrl_trgt_table: continue
            
            if ctrl_trgt_table[group_symbol] == 'control':
                group_0.append(expression)
            elif ctrl_trgt_table[group_symbol] == 'target':
                group_1.append(expression)

        GROUP_0[tx_id] = group_0
        GROUP_1[tx_id] = group_1

    for i in range(len(transcripts_data)):
        transcripts_data[i]['group0_exp'] = GROUP_0.get(transcripts_data[i]['id'], [])
        transcripts_data[i]['group1_exp'] = GROUP_1.get(transcripts_data[i]['id'], [])

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
                if len(col) >= 3:
                    tx_id = col[0]
                    logFC = float(col[1])
                    qval = float(col[2])
                    det_data[tx_id] = [logFC, qval]

    for i in range(len(transcripts_data)):
        if transcripts_data[i]['id'] in det_data:
            transcripts_data[i]['det'] = det_data[transcripts_data[i]['id']]
        else:
            transcripts_data[i]['det'] = [0, 1.0] # ダミー

    return transcripts_data


def prepare_de_ax1(transcripts_data, ax1, gene_name, tx_colors, colors_meta, tss_mode):
    # AX1: 各アイソフォームを視覚化
    y_positions = []
    if not transcripts_data: return ax1
    
    MIN = min(transcripts_data, key=lambda x: x['start'])['start']
    MAX = max(transcripts_data, key=lambda x: x['end'])['end']
    transcripts_data = transcripts_data[::-1]
    all_tss = []
    
    # 色情報の取得
    gcol = colors_meta.get("global", {})
    line_color = gcol.get("default_line", "gray")
    text_color = gcol.get("default_text", "black")
    default_tx_color = gcol.get("default_tx", "#B3C8CF")

    exon_override = colors_meta.get("exon", {}).get("color", None)
    cds_override = colors_meta.get("cds", {}).get("color", None)

    for i, transcript_data in enumerate(transcripts_data):
        y_positions.append(i)
        start = transcript_data['start']
        end = transcript_data['end']
        tx_name = transcript_data['id']
        strand = transcript_data['strand']

        # 1. 背骨 (イントロン)
        arrow_direction = "right" if strand == 1 else "left"
        y_pos = i
        
        # 背骨の線
        ax1.annotate('', xy=(end, y_pos), xytext=(start, y_pos), arrowprops=dict(arrowstyle="-", color=line_color, lw=1))
        
        # 矢印マーカー
        interval = max(1, (MAX-MIN)//50)
        x_positions = np.arange(start, end, interval)
        x_positions = x_positions[1:-1] # 両端を除く
        
        marker_shape = ">" if arrow_direction == "right" else "<"
        for x in x_positions:
            ax1.scatter(x, y_pos, marker=marker_shape, color=line_color, s=10)

        # TSS位置保存
        if strand == 1:
            all_tss.append(start)
        else:
            all_tss.append(end)

        # 色決定ロジック
        current_exon_color = exon_override or tx_colors.get(tx_name, default_tx_color)
        current_cds_color = cds_override or tx_colors.get(tx_name, default_tx_color) # FIX: tx_colorsを参照するように修正

        # 2. Exon を描画
        for exon in transcript_data['exons']:
            exon_start = exon[0]
            exon_end = exon[1]
            ax1.add_patch(patches.Rectangle((exon_start, i - 0.1), exon_end - exon_start, 0.2, color=current_exon_color))

        # 3. CDS を描画
        for cds in transcript_data['cds']:
            cds_start = cds[0]
            cds_end = cds[1]
            ax1.add_patch(patches.Rectangle((cds_start, i - 0.2), cds_end - cds_start, 0.4, color=current_cds_color))

    if tss_mode:
        unique_tss = sorted(set(all_tss))
        ymin = -0.3
        ymax = len(transcripts_data) - 0.8
        for x in unique_tss:
            ax1.vlines(x, ymin=ymin, ymax=ymax, colors=line_color, linestyles=':', linewidth=0.8, alpha=0.4)

    space = int((MAX-MIN)/20)
    ax1.set_xlim(MIN - space, MAX + space)
    ax1.set_ylim(-0.5, len(transcripts_data) - 0.5)
    
    chrom_label = transcripts_data[0]["seq_region_name"]
    chrom_norm = chrom_label.lower()
    if chrom_norm.startswith("chr"):
        chrom_norm = chrom_label[3:]
    else:
        chrom_norm = chrom_label

    ax1.set_xlabel(f'Chr{chrom_norm}', color=text_color)
    ax1.set_title(gene_name, color=text_color)
    ax1.set_yticks(y_positions)
    transcripts = [i['id'] for i in transcripts_data]
    ax1.set_yticklabels(transcripts, color=text_color)

    for label in ax1.get_xticklabels():
        label.set_rotation(15)
        label.set_ha('right')
        label.set_color(text_color)
        
    return ax1


def prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff):
    if ci == None:
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    else:
        ax1.set_xticks([x_min_eff, x_max_eff])
        ax1.set_xticklabels([str(x_min), str(x_max)])
    return ax1


def prepare_de_ax2(transcripts_data, ax2, config_meta, colors_meta, outliers, group0_label, group1_label):
    # AX2: 箱ひげ図
    transcripts_data = transcripts_data[::-1]
    transcripts = [i['id'] for i in transcripts_data]
    labels = []
    data = []
    
    # 配列順序: Group1(Target), Group0(Control)
    for i in range(len(transcripts_data)):
        tx_id = transcripts_data[i]['id']
        data.append(transcripts_data[i]['group1_exp']) # Index偶数 (0, 2...)
        data.append(transcripts_data[i]['group0_exp']) # Index奇数 (1, 3...)
        labels.append(f'{tx_id}_1')
        labels.append(f'{tx_id}_0')

    positions = []
    for i, name in enumerate(transcripts):
         positions.extend([i-0.2, i+0.2])

    boxplot = ax2.boxplot(data, positions=positions, vert=False, patch_artist=True, widths=0.2, showfliers=outliers)

    # 色を設定 (colors_metaから取得)
    c_group0 = colors_meta['de'].get('box_group0', '#4daf4a')
    c_group1 = colors_meta['de'].get('box_group1', '#e41a1c')
    
    box_colors = [c_group1 if idx % 2 == 0 else c_group0 for idx in range(len(data))]
    
    for patch, color in zip(boxplot['boxes'], box_colors):
        patch.set_facecolor(color)

    for median in boxplot['medians']:
        median.set_color('black')

    ax2.set_ylim(-0.5, len(transcripts_data) - 0.5)
    ax2_xlabel = "Expression level"
    ax2.set_yticklabels([])
    ax2.tick_params(axis='y', length=0)
    ax2.set_xlabel(ax2_xlabel)

    legend_handles = [
        patches.Patch(facecolor=c_group1, edgecolor='black', label=group1_label),
        patches.Patch(facecolor=c_group0, edgecolor='black', label=group0_label),
    ]
    ax2.legend(handles=legend_handles, frameon=False)
    return ax2


def prepare_de_ax3(transcripts_data, qval_threshold, ax3):
    # AX3: ヒートマップ (q値)
    def color_mapping(nlogqval):
        if nlogqval > 20:
            return 1
        else:
            return 0.3 + nlogqval*0.035
    
    DET = list()
    for i in range(len(transcripts_data)):
        DET.append([transcripts_data[i]['id']] + transcripts_data[i]['det'])

    warm_cmap = plt.cm.Reds
    cold_cmap = plt.cm.Blues

    for i, val in enumerate(DET):
        logFC = DET[i][1]
        qval = DET[i][2]
        
        if float(qval) < float(qval_threshold):
            if logFC > 0:
                color = warm_cmap(color_mapping(-np.log10(qval)))
            else:
                color = cold_cmap(color_mapping(-np.log10(qval)))
        else:
                color = (0.5, 0.5, 0.5, 0.5)

        rect = patches.Rectangle((0, len(DET) - i - 1), 0.1, 1, facecolor=color, edgecolor='black', linewidth=2)
        ax3.add_patch(rect)

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
        # orderにあるIDを優先して追加
        for tx in order:
            for i in transcripts_data:
                if i['id'] == tx:
                    tmp.append(i)
                    break
        # orderになかったものも残っていれば追加する場合（オプション）
        return tmp


def plot_isoespy_de(
        transcripts_data,
        config_meta,
        gene_name,
        tx_colors,
        colors_meta,
        ci,
        x_min,
        x_max,
        x_min_eff,
        x_max_eff,
        outliers,
        group0_label,
        group1_label,
        tss_mode,
        output_file=None,
        dpi=300,
        q_thresh=0.05):
    # 並び替え
    transcripts_data = reorder(transcripts_data, config_meta)
    
    if not transcripts_data:
        print("No transcripts found for plotting.")
        return

    # プロットする場所を作成
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(20, 5), gridspec_kw={'width_ratios': [10, 6, 2]})

    # prep for ax1 (colors_meta, tx_colorsを渡す)
    ax1 = prepare_de_ax1(transcripts_data, ax1, gene_name, tx_colors, colors_meta, tss_mode)
    ax1 = prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff)

    # prep for ax2 (colors_metaを渡す)
    ax2 = prepare_de_ax2(transcripts_data, ax2, config_meta, colors_meta, outliers, group0_label, group1_label)

    # prep for ax3
    ax3 = prepare_de_ax3(transcripts_data, q_thresh, ax3)

    plt.subplots_adjust(wspace=0.02)

    if output_file:
        try:
            fig.savefig(
                output_file, 
                dpi=dpi, 
                bbox_inches='tight'
            )
            print(f"Figure successfully saved to {output_file} with DPI={dpi}")
        except Exception as e:
            print(f"Error saving figure to file: {e}")
            plt.show() # fallback
    else:
        plt.show()


def moved_data_for_exons_cds(main_data, ci):
    """ intron compression 処理 """
    # preparation
    model = {}
    for tx_data in main_data:
        tx_name = tx_data['id']
        model[tx_name+"_exons"] = tx_data['exons']
        model[tx_name+"_cds"] = tx_data['cds']

    model_compressed = intronCompression(model, ci)

    startend_d = {}
    for tx_data in main_data:
        tx_name = tx_data['id']
        if tx_name+"_exons" in model_compressed:
            exons = model_compressed[tx_name+"_exons"]
            if exons:
                start = exons[0][0]
                end = exons[-1][1]
                startend_d[tx_name] = {"start": start, "end": end}
                tx_data['exons'] = exons
                tx_data['cds'] = model_compressed[tx_name+"_cds"]
                tx_data['start'] = start
                tx_data['end'] = end

    return main_data


def isoespy_de(gene, gtf_data, expression_data, det_data, meta_data, ci, outliers, group0_label, group1_label, tss_mode, output_file=None, dpi=300, q_thresh=0.05):
    # metadata
    # parse_metadata が返す値を更新（style_meta -> colors_meta）
    sample_meta, config_meta, ctrl_trgt_table, gtf_meta, colors_meta, query, feature_meta = parse_metadata(meta_data, gene)

    # raw isoform model (colors_meta, tx_colorsを処理)
    tx_colors = colors_meta["transcripts"]
    model_exon, model_cds, tx_colors = get_isoform_model(gtf_data, gtf_meta, tx_colors, query, colors_meta)

    # formatted isoform model
    transcripts_main_data = list()
    transcripts_main_data = formatting_isoform_model(transcripts_main_data, model_exon, annot="exon")
    transcripts_main_data = formatting_isoform_model(transcripts_main_data, model_cds, annot="cds")
    
    if not transcripts_main_data:
        print(f"No transcripts found for gene: {gene}")
        return

    x_min = min(transcripts_main_data, key=lambda x: x['start'])['start']
    x_max = max(transcripts_main_data, key=lambda x: x['end'])['end']

    if ci != None:
        transcripts_main_data = moved_data_for_exons_cds(transcripts_main_data, ci)
    
    x_min_eff = min(transcripts_main_data, key=lambda x: x['start'])['start']
    x_max_eff = max(transcripts_main_data, key=lambda x: x['end'])['end']

    transcripts_main_data = get_expression_data(expression_data, transcripts_main_data, ctrl_trgt_table, sample_meta)
    transcripts_main_data = get_det_data(det_data, transcripts_main_data)

    plot_isoespy_de(
        transcripts_main_data,
        config_meta,
        gene,
        tx_colors,   # 追加
        colors_meta, # style_meta -> colors_meta
        ci,
        x_min,
        x_max,
        x_min_eff,
        x_max_eff,
        outliers,
        group0_label,
        group1_label,
        tss_mode,
        output_file,
        dpi,
        q_thresh
    )


def process_ci(ci):
    if ci is None: return None
    if isinstance(ci, str):
        try:
            ci = float(ci)
        except ValueError:
            warnings.warn(f"Warning: ci should be a float/int, but received string '{ci}'")
            return None
    return ci


def main(args=None):
    parser = argparse.ArgumentParser(description='isoespy_de()')
    parser.add_argument('-gene', '--gene_name', required=True, type=str, default=None, help='Gene name')
    parser.add_argument('-gtf', '--gtf_data', required=True, type=str, default=None, help='GTF file')
    parser.add_argument('-exp', '--expression_data', required=True, type=str, default=None, help='Expression data')
    parser.add_argument('-det', '--det_data', required=True, type=str, default=None, help='Differential expression data')
    parser.add_argument('-meta', '--meta_data', required=True, type=str, default=None, help='metadata')
    parser.add_argument('-ci', '--compress_introns', default=None, help='intron compression parameter')
    parser.add_argument('--show_outliers', action="store_true", help="Show outliers in the boxplot (default: hide)")
    parser.add_argument('--group0_label', default='0', help='Label of control (eg. Nontumor)')
    parser.add_argument('--group1_label', default='1', help='Label of target (eg. HCC)')
    parser.add_argument('-tss', '--tss_line', action="store_true", help="Show TSS support lines")
    parser.add_argument('-qval', '--q_thresh', type=float, default=0.05, help='q-value/FDR threshold for coloring significant isoforms')
    parser.add_argument("-o", "--output_file", type=str, default=None, help="Save figure to file (e.g., plot.pdf or plot.png)")
    parser.add_argument("--dpi", type=int, default=300, help="Resolution for raster image formats (e.g., PNG) in DPI")

    args = parser.parse_args()
    gene = args.gene_name
    gtf_data = args.gtf_data
    expression_data = args.expression_data
    det_data = args.det_data
    meta_data = args.meta_data
    ci = process_ci(args.compress_introns)
    outliers = args.show_outliers
    group0_label = args.group0_label
    group1_label = args.group1_label
    tss_mode = args.tss_line
    q_thresh = args.q_thresh
    output_file = args.output_file
    dpi = args.dpi

    isoespy_de(gene, gtf_data, expression_data, det_data, meta_data, ci, outliers, group0_label, group1_label, tss_mode, output_file, dpi, q_thresh)


if __name__ == '__main__':
    main()
