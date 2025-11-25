import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import warnings
import argparse
import matplotlib.ticker as ticker
import math
import re
import seaborn as sns  # 追加
import matplotlib.colors as mcolors # 追加

try:
    from intronCompression import intronCompression
except ModuleNotFoundError:
    try:
        from isoespy.intronCompression import intronCompression
    except ModuleNotFoundError:
        def intronCompression(model, ci):
            return model

import matplotlib.pyplot as plt
plt.rcParams.update({
    'font.size': 8,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
})

# -----------------------------
# metadata parser (Refactored to match isoespy_ff/de)
# -----------------------------
def parse_metadata(meta_data, gene):
    """
    メタデータファイルをパースして色設定(colors_meta)などを含めて返す。
    isoespy_ff/de と同じ階層構造を採用。
    """
    with open(meta_data) as f:
        lines = f.readlines()

    current_section = None
    gtf_meta = {}
    config_meta = {}
    query = {"gene": [None, None], "tx": [None, None]}
    
    # Hierarchical colors meta
    colors_meta = {
        "global": {
            "palette": None,
            "default_line": "gray",
            "default_text": "black",
            "default_tx": "#B3C8CF",
        },
        "transcripts": {},
        "exon": {"color": None},
        "cds": {"color": None},
        # DIU specific settings
        "diu": {
            "bar_group0": "#fdd9b5", # Default control color
            "bar_group1": "#d95f02"  # Default target color
        }
    }
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith("!"): 
             continue

        if line.startswith("[") and line.endswith("]"):
            current_section = line[1:-1]
            continue

        # パラメータ行 (#key=value または key=value)
        content = line.lstrip("#").strip() if line.startswith("#") else line

        if current_section == "config":
            if "=" in content:
                key, value = content.split("=", 1)
                key = key.strip()
                value = value.strip()
                if key == 'order':
                    config_meta['order'] = [i.strip() for i in value.split(',')]
                elif key == 'colors':
                    # Legacy style support: colors=<color>:tx1,tx2...
                    if ":" in value:
                        val1, val2 = value.split(":", 1)
                        val1 = val1.strip()
                        val2 = [i.strip() for i in val2.split(",")]
                        for tx in val2:
                             colors_meta['transcripts'][tx.strip()] = val1
                else:
                    config_meta[key] = value

        elif current_section == "gtf":
            if "=" in content:
                key, value = content.split("=", 1)
                gtf_meta[key.strip()] = value.strip()

        elif current_section == "query":
            if "=" in content:
                key, value = content.split("=", 1)
                key = key.strip()
                value = value.strip()
                if key == 'gene':
                    query["gene"] = [value, gene]
                elif key == 'tx':
                    tx_list = [i.strip() for i in value.split(",")]
                    query["tx"][1] = set(tx_list)
        
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
                colors_meta['transcripts'][key.strip()] = value.strip()
        
        # ---------------- colors.exon / cds ----------------
        elif current_section == "colors.exon":
            if "=" in content:
                key, value = content.split("=", 1)
                if key.strip() == "color":
                    colors_meta['exon']['color'] = value.strip()

        elif current_section == "colors.cds":
            if "=" in content:
                key, value = content.split("=", 1)
                if key.strip() == "color":
                    colors_meta['cds']['color'] = value.strip()
        
        # ---------------- colors.diu (Specific) ----------------
        elif current_section == "colors.diu":
            if "=" in content:
                key, value = content.split("=", 1)
                key = key.strip()
                if key in colors_meta['diu']:
                    colors_meta['diu'][key] = value.strip()

    query['tx'][0] = gtf_meta.get('transcript_id')
    if query['tx'][1] == set() or query['tx'][1] == None:
        query['tx'][1] = None

    return config_meta, gtf_meta, colors_meta, query


def get_attr_candidate(attr_dict, keys):
    for k in keys:
        if k in attr_dict:
            return attr_dict[k]
    return None


def get_isoform_model(gtf_file, gtf_meta, tx_colors, query, colors_meta=None):
    """
    GTFからアイソフォーム構造を取得。色の自動補完も行う。
    """
    transcripts = {}
    transcripts_CDS = {}

    target_gene = query["gene"][1]
    target_tx_filter = query["tx"][1]
    tx_id_key_from_meta = gtf_meta.get('transcript_id', 'transcript_id')

    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            col = line.strip().split('\t')
            if len(col) < 9:
                continue

            chrom, _, ftype, start, end, _, strand, _, attr = col
            start = int(start)
            end = int(end)
            strand = 1 if strand == "+" else -1

            attr_dict = {}
            for match in re.finditer(r'(\S+)\s+"([^"]+)"', attr):
                attr_dict[match.group(1)] = match.group(2)
            if not attr_dict:
                 for raw in attr.strip().split(";"):
                    parts = raw.strip().split("=") 
                    if len(parts) == 2:
                        attr_dict[parts[0].strip()] = parts[1].strip().strip('"')

            line_gene = get_attr_candidate(attr_dict, [gtf_meta.get('gene_name', 'gene_name'), 'gene_name', 'gene_symbol', 'gene_id', 'Name'])
            line_tx = get_attr_candidate(attr_dict, [tx_id_key_from_meta, 'transcript_id', 'transcript', 'tx_id', 'ID'])

            if line_gene != target_gene:
                continue

            if target_tx_filter is not None:
                if line_tx not in target_tx_filter:
                    continue

            if ftype == "exon":
                if line_tx not in transcripts:
                    transcripts[line_tx] = [[], strand, chrom]
                transcripts[line_tx][0].append((start, end))

            if ftype == "CDS":
                if line_tx not in transcripts_CDS:
                    transcripts_CDS[line_tx] = [[], strand, chrom]
                transcripts_CDS[line_tx][0].append((start, end))

    for tx in transcripts:
        transcripts[tx][0] = sorted(transcripts[tx][0], key=lambda x: x[0])
    for tx in transcripts_CDS:
        transcripts_CDS[tx][0] = sorted(transcripts_CDS[tx][0], key=lambda x: x[0])

    if len(transcripts) == 0:
        print("[WARN] No transcripts found for gene:", target_gene)

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
            end   = max([i[1] for i in exons[0]])
            transcripts_data.append({
                'id': transcript_id,
                'exons': exons[0],
                'strand': exons[1],
                'seq_region_name': exons[2],
                'start': start,
                'end': end,
            })
    elif annot == "cds":
        for i in range(len(transcripts_data)):
            tx_id = transcripts_data[i]['id']
            if tx_id in transcripts:
                transcripts_data[i]['cds'] = transcripts[tx_id][0]
            else:
                transcripts_data[i]['cds'] = []
    return transcripts_data


def moved_data_for_exons_cds(main_data, ci):
    model = {}
    for tx_data in main_data:
        tx = tx_data['id']
        model[tx+"_exons"] = tx_data['exons']
        model[tx+"_cds"]   = tx_data['cds']

    model_compressed = intronCompression(model, ci)

    startend_d = {}
    for tx_data in main_data:
        tx = tx_data['id']
        if tx+"_exons" in model_compressed and model_compressed[tx+"_exons"]:
            startend_d[tx] = {
                "start": model_compressed[tx+"_exons"][0][0],
                "end":   model_compressed[tx+"_exons"][-1][1],
            }
        else:
            startend_d[tx] = {"start": tx_data['start'], "end": tx_data['end']}

    for tx_data in main_data:
        tx = tx_data['id']
        tx_data['exons'] = model_compressed[tx+"_exons"]
        tx_data['cds']   = model_compressed[tx+"_cds"]
        tx_data['start'] = startend_d[tx]["start"]
        tx_data['end']   = startend_d[tx]["end"]

    return main_data


def prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff):
    if ci is None:
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    else:
        ax1.set_xticks([x_min_eff, x_max_eff])
        ax1.set_xticklabels([str(x_min), str(x_max)])
    return ax1


def prepare_ax1_isoform_structure(transcripts_data, ax1, gene_name, tx_colors, colors_meta, tss_mode):
    # AX1: 構造図 (色設定を反映)
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

    y_positions = []
    for i, transcript_data in enumerate(transcripts_data):
        y_positions.append(i)
        start = transcript_data['start']
        end   = transcript_data['end']
        tx_id = transcript_data['id']
        strand = transcript_data['strand']
        
        arrow_direction = "right" if strand == 1 else "left"
        y_pos = i
        
        # 背骨
        ax1.annotate('', xy=(end, y_pos), xytext=(start, y_pos),
                     arrowprops=dict(arrowstyle="-", color=line_color, lw=1))
        
        interval = max(1, (MAX - MIN)//50)
        x_positions = np.arange(start, end, interval)[1:-1]
        
        if arrow_direction == "right":
            for x in x_positions: ax1.scatter(x, y_pos, marker=">", color=line_color, s=10)
            all_tss.append(start)
        else:
            for x in x_positions: ax1.scatter(x, y_pos, marker="<", color=line_color, s=10)
            all_tss.append(end)

        # 色決定
        current_exon_color = exon_override or tx_colors.get(tx_id, default_tx_color)
        current_cds_color = cds_override or tx_colors.get(tx_id, default_tx_color) # FIX: tx_colorsを参照するように修正

        # exon
        for exon_start, exon_end in transcript_data['exons']:
            ax1.add_patch(patches.Rectangle((exon_start, i - 0.1), exon_end - exon_start, 0.2, color=current_exon_color))

        # CDS
        for cds_start, cds_end in transcript_data['cds']:
            ax1.add_patch(patches.Rectangle((cds_start, i - 0.2), cds_end - cds_start, 0.4, color=current_cds_color))

    space = int((MAX-MIN)/20)
    ax1.set_xlim(MIN - space, MAX + space)
    ax1.set_ylim(-0.5, len(transcripts_data) - 0.5)

    if tss_mode:
        unique_tss = sorted(set(all_tss))
        ymin = -0.3
        ymax = len(transcripts_data) - 0.8
        for x in unique_tss:
            ax1.vlines(x, ymin=ymin, ymax=ymax, colors=line_color, linestyles=':', linewidth=0.8, alpha=0.4)

    # Chromosome label formatting
    chrom_label = transcripts_data[0]["seq_region_name"]
    chrom_norm = chrom_label.lower()
    if chrom_norm.startswith("chr"):
        chrom_norm = chrom_label[3:]
    else:
        chrom_norm = chrom_label

    ax1.set_xlabel(f'Chr{chrom_norm}', color=text_color)
    ax1.set_title(gene_name, color=text_color)

    ax1.set_yticks(y_positions)
    ax1.set_yticklabels([t['id'] for t in transcripts_data], color=text_color)

    for label in ax1.get_xticklabels():
        label.set_rotation(15)
        label.set_ha('right')
        label.set_color(text_color)

    return ax1


def get_diu_data(diu_file, transcripts_data, gene_symbol,
                 col_geneID="geneID",
                 col_tx="feature_id",
                 col_prop_ctrl="prop_N",
                 col_prop_case="prop_T",
                 col_p="transcript"):
    
    df = pd.read_csv(diu_file, sep="\t")
    # geneID列でフィルタリング
    if col_geneID in df.columns:
        df_gene = df[df[col_geneID] == gene_symbol].copy()
    else:
        # カラム名が見つからない場合のフォールバック（またはエラー回避）
        df_gene = df
    
    if df_gene.empty:
         # データがない場合はNaNで埋める
         for t in transcripts_data:
            t['group0_prop'] = np.nan
            t['group1_prop'] = np.nan
            t['delta_prop']  = np.nan
            t['pval_diu']    = np.nan
         return transcripts_data

    df_gene_indexed = df_gene.set_index(col_tx)

    for t in transcripts_data:
        tx = t['id']
        if tx in df_gene_indexed.index:
            row = df_gene_indexed.loc[tx]
            try:
                ctrl_prop = float(row[col_prop_ctrl])
                case_prop = float(row[col_prop_case])
                delta = case_prop - ctrl_prop
                pval = float(row[col_p])
            except (KeyError, ValueError):
                ctrl_prop, case_prop, delta, pval = np.nan, np.nan, np.nan, np.nan
        else:
            ctrl_prop, case_prop, delta, pval = np.nan, np.nan, np.nan, np.nan

        t['group0_prop'] = ctrl_prop
        t['group1_prop'] = case_prop
        t['delta_prop']  = delta
        t['pval_diu']    = pval

    return transcripts_data


def reorder(transcripts_data, order_list):
    if not order_list:
        return transcripts_data
    out = []
    # orderにあるものを先に
    for tx in order_list:
        for t in transcripts_data:
            if t['id'] == tx:
                out.append(t)
                break
    # orderにないものも表示したい場合は以下を追加（今回はorderのみ抽出する挙動とする）
    return out


def prepare_diu_ax2(transcripts_data, ax2, colors_meta,
                    group0_label="Nontumor",
                    group1_label="HCC"):
    """
    AX2: Mean Usage Bar Plot (色指定反映)
    """
    tdata = transcripts_data[::-1]
    y_pos = np.arange(len(tdata))

    group0_vals = [t['group0_prop'] for t in tdata]
    group1_vals = [t['group1_prop'] for t in tdata]

    # メタデータから色を取得 (colors_meta['diu']から)
    c_group0 = colors_meta['diu'].get('bar_group0', '#fdd9b5') # control
    c_group1 = colors_meta['diu'].get('bar_group1', '#d95f02') # target

    bar_h = 0.35
    ax2.barh(y_pos - bar_h/2, group1_vals, height=bar_h, label=group1_label, color=c_group1)
    ax2.barh(y_pos + bar_h/2, group0_vals, height=bar_h, label=group0_label, color=c_group0)

    # usage max
    max_usage = max(np.nanmax(group1_vals), np.nanmax(group0_vals))
    if not np.isfinite(max_usage) or max_usage <= 0:
        max_usage = 1.0

    ax2.set_xlim(0, max_usage * 1.05)
    ax2.set_ylim(-0.5, len(tdata) - 0.5)
    ax2.set_yticks([])
    ax2.tick_params(axis='y', length=0)

    from matplotlib.ticker import FuncFormatter
    ax2.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '0' if x == 0 else f'{x:.2f}'))

    ax2.set_xlabel("Mean usage")
    ax2.legend(fontsize=8, frameon=False)
    return ax2


def prepare_diu_ax3(transcripts_data, ax3,
                    p_thresh=0.05,
                    group0_label="Nontumor",
                    group1_label="HCC"):
    tdata = transcripts_data[::-1]
    y_pos = np.arange(len(tdata))

    delta_vals = np.array([t['delta_prop'] for t in tdata], dtype=float)
    pvals      = np.array([t['pval_diu']   for t in tdata], dtype=float)

    sig_mask = (pvals < p_thresh) & (~np.isnan(pvals)) & (~np.isnan(delta_vals))

    with np.errstate(divide='ignore'):
        logp = -np.log10(pvals)

    if np.any(sig_mask):
        max_logp_sig = np.nanmax(logp[sig_mask])
        if (not np.isfinite(max_logp_sig)) or max_logp_sig <= 0:
            max_logp_sig = 1.0
    else:
        max_logp_sig = 1.0

    colors = []
    for du, pv, lp, sig in zip(delta_vals, pvals, logp, sig_mask):
        if (not sig) or (not np.isfinite(du)):
            colors.append((0.7,0.7,0.7)) # Gray
        else:
            intensity = lp / max_logp_sig
            intensity = max(0.0, min(1.0, intensity))
            if du >= 0:
                # Red scale
                r, g, b = 1.0, 0.8 - 0.5*intensity, 0.8 - 0.5*intensity
                colors.append((r,g,b))
            else:
                # Blue scale
                r, g, b = 0.8 - 0.5*intensity, 0.8 - 0.5*intensity, 1.0
                colors.append((r,g,b))

    max_abs = np.nanmax(np.abs(delta_vals)) if len(delta_vals) > 0 else 0.1
    if (not np.isfinite(max_abs)) or max_abs == 0: max_abs = 0.1
    pad = max_abs * 0.2
    xlim_delta = (-max_abs - pad, max_abs + pad)

    ax3.barh(y_pos, delta_vals, height=0.6, color=colors, edgecolor="none")
    ax3.axvline(0, color="k", linewidth=0.8)
    ax3.set_xlim(xlim_delta)
    ax3.set_ylim(-0.5, len(tdata) - 0.5)

    from matplotlib.ticker import FuncFormatter
    ax3.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '0' if x == 0 else f'{x:.2f}'))
    ax3.set_yticks([])
    ax3.tick_params(axis='y', length=0)
    ax3.set_xlabel(f"Δusage ({group1_label} - {group0_label})")

    return ax3


def process_ci(ci):
    if ci is None: return None
    try:
        return float(ci)
    except (ValueError, TypeError):
        warnings.warn(f"Warning: ci should be numeric, got {ci}. intron compression disabled.")
        return None


def isoespy_diu(gene_symbol,
                gtf_data,
                diu_tsv,
                meta_data,
                ci=None,
                group0_label="0",
                group1_label="1",
                p_thresh=0.05,
                tss_mode=False,
                output_file=None,
                dpi=300):
    
    # 1. メタ情報読込 (colors_metaを取得)
    config_meta, gtf_meta, colors_meta, query = parse_metadata(meta_data, gene_symbol)

    # 2. GTF読込 (tx_colorsとcolors_metaを渡す)
    tx_colors = colors_meta['transcripts']
    transcripts_exon, transcripts_cds, tx_colors = get_isoform_model(gtf_data, gtf_meta, tx_colors, query, colors_meta)

    # 3. フォーマット
    transcripts_main = []
    transcripts_main = formatting_isoform_model(transcripts_main, transcripts_exon, annot="exon")
    transcripts_main = formatting_isoform_model(transcripts_main, transcripts_cds,  annot="cds")

    if not transcripts_main:
        print(f"No transcripts found for gene: {gene_symbol}")
        return

    x_min = min(transcripts_main, key=lambda x: x['start'])['start']
    x_max = max(transcripts_main, key=lambda x: x['end'])['end']

    # 4. イントロン圧縮
    ci_val = process_ci(ci)
    if ci_val is not None:
        transcripts_main = moved_data_for_exons_cds(transcripts_main, ci_val)

    x_min_eff = min(transcripts_main, key=lambda x: x['start'])['start']
    x_max_eff = max(transcripts_main, key=lambda x: x['end'])['end']

    # 5. DIUデータを付加
    transcripts_main = get_diu_data(diu_tsv, transcripts_main, gene_symbol)

    # 6. 並び順
    transcripts_main = reorder(transcripts_main, config_meta.get('order'))
    
    if not transcripts_main:
        print("No transcripts left after filtering.")
        return

    # 7. 描画
    fig, (ax1, ax2, ax3) = plt.subplots(
        ncols=3, figsize=(20,5), gridspec_kw={'width_ratios':[10,6,3]}
    )

    # ax1: 構造 (tx_colors, colors_metaを渡す)
    ax1 = prepare_ax1_isoform_structure(transcripts_main, ax1, gene_symbol, tx_colors, colors_meta, tss_mode)
    ax1 = prepare_ax1_xaxis(ax1, ci_val, x_min, x_max, x_min_eff, x_max_eff)

    # ax2: Usage (colors_metaを渡す)
    ax2 = prepare_diu_ax2(transcripts_main, ax2, colors_meta,
                          group0_label=group0_label,
                          group1_label=group1_label)

    # ax3: Delta Usage
    ax3 = prepare_diu_ax3(transcripts_main, ax3,
                          p_thresh=p_thresh,
                          group0_label=group0_label,
                          group1_label=group1_label)

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


def main():
    parser = argparse.ArgumentParser(description="isoespy_diu: plot isoform structure + DIU usage/Δusage panels.")
    parser.add_argument('-gene', '--gene_name', required=True, help="Gene name")
    parser.add_argument('-gtf', '--gtf_data', required=True, help="GTF file")
    parser.add_argument('-diu', "--diu_result", required=True, help="DIU results TSV")
    parser.add_argument('-meta', '--meta_data', required=True, help="metadata file")
    parser.add_argument('-ci', '--compress_introns', default=None, help="float. Apply intron compression with this factor")
    parser.add_argument('--group0_label', default='0', help="label for control/baseline group (prop_N)")
    parser.add_argument('--group1_label', default='1', help="label for target/case group (prop_T)")
    parser.add_argument('-pval', '--p_thresh', default=0.05, type=float, help="q-value/p-value cutoff for coloring Δusage")
    parser.add_argument('-tss', '--tss_line', action="store_true", help="Show TSS support lines")
    parser.add_argument("-o", "--output_file", type=str, default=None, help="Save figure to file (e.g., plot.pdf or plot.png)")
    parser.add_argument("--dpi", type=int, default=300, help="Resolution for raster image formats (e.g., PNG) in DPI")

    args = parser.parse_args()

    isoespy_diu(
        gene_symbol=args.gene_name,
        gtf_data=args.gtf_data,
        diu_tsv=args.diu_result,
        meta_data=args.meta_data,
        ci=args.compress_introns,
        group0_label=args.group0_label,
        group1_label=args.group1_label,
        p_thresh=args.p_thresh,
        tss_mode=args.tss_line,
        output_file=args.output_file,
        dpi=args.dpi
    )

if __name__ == "__main__":
    main()
