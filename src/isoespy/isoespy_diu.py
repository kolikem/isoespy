import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import warnings
import argparse
import matplotlib.ticker as ticker
import math
import re
import seaborn as sns
import matplotlib.colors as mcolors
import copy
from collections import defaultdict

try:
    from intronCompression import intronCompression
except ModuleNotFoundError:
    try:
        from isoespy.intronCompression import intronCompression
    except ModuleNotFoundError:
        def intronCompression(model, ci):
            return model

# Plot Style Settings
plt.rcParams.update({
    'font.size': 8,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans']
})


# Helper Functions
def get_attr_candidate(attr_dict, keys):
    """Try multiple keys to find a value in attr_dict."""
    for k in keys:
        if k in attr_dict:
            return attr_dict[k]
    return None

def process_ci(ci):
    """Parse the compress_introns parameter."""
    if ci is None: return None
    try:
        return float(ci)
    except (ValueError, TypeError):
        warnings.warn(f"Warning: ci should be numeric, got {ci}. Intron compression disabled.")
        return None

def darken_color(color, factor=0.4):
    return tuple(max(0, c * factor) for c in color)


# Metadata Parser
def parse_metadata(meta_data, gene):
    with open(meta_data) as f:
        lines = f.readlines()

    current_section = None
    gtf_meta = {}
    config_meta = {}
    query = {"gene": [None, None], "tx": [None, None]}
    feature_meta = dict()
    
    # Hierarchical colors meta
    colors_meta = {
        "global": {
            "palette": None,
            "default_line": "gray",
            "default_text": "black",
            "default_tx": "#999999",
            "default_cds": "#3b73b9",
        },
        "transcripts": {},
        "exon": {"color": None},
        "cds": {"color": None},
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

        content = line.lstrip("#").strip() if line.startswith("#") else line

        if current_section == "config":
            if "=" in content:
                key, value = content.split("=", 1)
                key = key.strip()
                value = value.strip()
                if key == 'order':
                    config_meta['order'] = [i.strip() for i in value.split(',')]
                elif key == 'colors':
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

        elif current_section == "features":
            if ":" in content:
                feature_class, feature_id = content.split(":", 1)
                feature_meta[feature_class.strip()] = feature_id.strip() if feature_id.strip() else None
            else:
                feature_meta[content] = None
        
        elif current_section == "colors":
            if "=" in content:
                key, value = content.split("=", 1)
                if key.strip() in colors_meta["global"]:
                    colors_meta["global"][key.strip()] = value.strip()

        elif current_section == "colors.transcripts":
            if "=" in content:
                key, value = content.split("=", 1)
                colors_meta['transcripts'][key.strip()] = value.strip()
        
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
        
        elif current_section == "colors.diu":
            if "=" in content:
                key, value = content.split("=", 1)
                key = key.strip()
                if key in colors_meta['diu']:
                    colors_meta['diu'][key] = value.strip()

    query['tx'][0] = gtf_meta.get('transcript_id')
    if query['tx'][1] == set() or query['tx'][1] == None:
        query['tx'][1] = None

    return config_meta, gtf_meta, colors_meta, query, feature_meta


def get_scale_unit(data_range):
    target_size = data_range * 0.15
    order_of_magnitude = 10 ** int(np.floor(np.log10(target_size)))
    multipliers = [1, 2, 5, 10]
    selected_size = order_of_magnitude
    for m in multipliers:
        if order_of_magnitude * m <= target_size:
            selected_size = order_of_magnitude * m
        else:
            break
    if selected_size >= 1000:
        label = f"{int(selected_size/1000)} kb"
    else:
        label = f"{int(selected_size)} bp"
    return selected_size, label


# ---------------------------------------------------------
# Model Builders
# ---------------------------------------------------------
def get_isoform_model(gtf_file, gtf_meta, tx_colors, query, colors_meta=None):
    transcripts = {}
    transcripts_CDS = {}
    target_gene = query["gene"][1]
    target_tx_filter = query["tx"][1]
    tx_id_key_from_meta = gtf_meta.get('transcript_id', 'transcript_id')

    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith('#'): continue
            col = line.strip().split('\t')
            if len(col) < 9: continue

            chrom, _, ftype, start, end, _, strand, _, attr = col
            start = int(start)
            end = int(end)
            strand = 1 if strand == "+" else -1

            attr_dict = {match.group(1): match.group(2) for match in re.finditer(r'(\S+)\s+"([^"]+)"', attr)}
            if not attr_dict:
                 for raw in attr.strip().split(";"):
                    parts = raw.strip().split("=") 
                    if len(parts) == 2:
                        attr_dict[parts[0].strip()] = parts[1].strip().strip('"')

            line_gene = get_attr_candidate(attr_dict, [gtf_meta.get('gene_name', 'gene_name'), 'gene_name', 'gene_symbol', 'gene_id', 'Name'])
            line_tx = get_attr_candidate(attr_dict, [tx_id_key_from_meta, 'transcript_id', 'transcript', 'tx_id', 'ID'])

            if line_gene != target_gene: continue
            if target_tx_filter is not None:
                if line_tx not in target_tx_filter: continue

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
            for tx in missing:
                tx_colors[tx] = default_tx
    
    return transcripts, transcripts_CDS, tx_colors


def get_feature_model(gtf_file, gtf_meta, query, feature_meta):
    ff_d = dict()
    target_gene = query["gene"][1]
    target_tx = query["tx"][1]
    gene_keys = [query["gene"][0], 'gene_name', 'gene_symbol', 'gene_id', 'Name']
    gene_keys = [k for k in gene_keys if k]

    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith("#"): continue
            fields = line.strip().split("\t")
            if len(fields) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            if feature in feature_meta:
                attr_dict = {match.group(1): match.group(2) for match in re.finditer(r'(\S+)\s+"([^"]+)"', attributes)}
                if not attr_dict:
                    for raw in attributes.strip().split(";"):
                        parts = raw.strip().split("=") 
                        if len(parts) == 2:
                            attr_dict[parts[0].strip()] = parts[1].strip().strip('"')

                transcript_id = get_attr_candidate(attr_dict, [gtf_meta.get("transcript_id"), "transcript_id", "transcript", "ID"])
                line_gene = get_attr_candidate(attr_dict, gene_keys)
                line_tx = get_attr_candidate(attr_dict, [query["tx"][0], "transcript_id", "transcript", "ID"])

                match_gene = (line_gene == target_gene)
                match_tx = (target_tx is None) or (line_tx in target_tx)
                if not match_gene and line_tx and target_tx and (line_tx in target_tx):
                    match_gene = True

                if not match_gene: continue
                if not match_tx: continue

                if transcript_id not in ff_d:
                    ff_d[transcript_id] = {feat: {} for feat in feature_meta}

                if feature_meta[feature] is None:
                    if ff_d[transcript_id][feature] == {}:
                        ff_d[transcript_id][feature]["NONAME"] = []
                    ff_d[transcript_id][feature]["NONAME"].append((int(start), int(end)))
                else:
                    feat_id = attr_dict.get(feature_meta[feature])
                    if feat_id:
                        if feat_id not in ff_d[transcript_id][feature]:
                            ff_d[transcript_id][feature][feat_id] = []
                        ff_d[transcript_id][feature][feat_id].append((int(start), int(end)))

    for tx in ff_d:
        for feat in feature_meta:
            for ind in ff_d[tx][feat]:
                ff_d[tx][feat][ind] = sorted(ff_d[tx][feat][ind], key=lambda x: x[0])
    return ff_d


def coordinate_adjustment(feat_names, Y_LOW, Y_UPP):
    out_d = dict()
    if feat_names == []: return out_d
    width = (Y_UPP - Y_LOW) / len(feat_names)
    for i in range(len(feat_names)):
        out_d[feat_names[i]] = (Y_UPP - width * (i + 1), Y_UPP - width * i)
    return out_d


def overlaped_feature_indivisuals(DICT):
    def overlaps(r1, r2):
        return r1[0] <= r2[1] and r2[0] <= r1[1]

    graph = defaultdict(set)
    keys = list(DICT.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
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


def moved_data_for_features(model_features, main_data, ci):
    model = {}
    for tx_data in main_data:
        tx_name = tx_data["id"]
        model[tx_name + "_exons"] = tx_data["exons"]
        model[tx_name + "_cds"] = tx_data["cds"]

    for tx_name in model_features:
        for feat_name in model_features[tx_name]:
            for id_name in model_features[tx_name][feat_name]:
                model[tx_name + "_" + feat_name + "_" + id_name] = model_features[tx_name][feat_name][id_name]

    model_compressed = intronCompression(model, ci)

    for tx_name in model_features:
        for feat_name in model_features[tx_name]:
            for id_name in model_features[tx_name][feat_name]:
                key = tx_name + "_" + feat_name + "_" + id_name
                if key in model_compressed:
                    model_features[tx_name][feat_name][id_name] = model_compressed[key]
    return model_features


# ---------------------------------------------------------
# Plotting Functions
# ---------------------------------------------------------

def prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff):
    if ci is None:
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    else:
        ax1.set_xticks([x_min_eff, x_max_eff])
        ax1.set_xticklabels([str(x_min), str(x_max)])
    return ax1


def prepare_ax1_isoform_structure(transcripts_data, ax1, gene_name, tx_colors, colors_meta, tss_mode, model_features, feature_meta, ci=None, hide_scale=False):
    # AX1: Structure
    if not transcripts_data: return ax1
    
    MIN = min(transcripts_data, key=lambda x: x['start'])['start']
    MAX = max(transcripts_data, key=lambda x: x['end'])['end']
    transcripts_data = transcripts_data[::-1]
    all_tss = []
    
    if model_features == {}:
        feat_names = []
        feat_colors = []
    else:
        all_found_keys = set()
        for tx in model_features:
            all_found_keys.update(model_features[tx].keys())
        feat_names = [k for k in feature_meta.keys() if k in all_found_keys]
        feat_colors = sns.color_palette("husl", len(feat_names))

    gcol = colors_meta.get("global", {})
    line_color = gcol.get("default_line", "gray")
    text_color = gcol.get("default_text", "black")
    default_tx_color = gcol.get("default_tx", "#999999")
    default_cds_color = gcol.get("default_cds", "#3b73b9")
    exon_override = colors_meta.get("exon", {}).get("color", None)
    cds_override = colors_meta.get("cds", {}).get("color", None)
    user_tx_colors = colors_meta.get("user_transcripts", {})

    y_positions = []
    for i, transcript_data in enumerate(transcripts_data):
        y_positions.append(i)
        start = transcript_data['start']
        end   = transcript_data['end']
        tx_id = transcript_data['id']
        strand = transcript_data['strand']
        
        arrow_direction = "right" if strand == 1 else "left"
        y_pos = i
        
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

        # Exon Color
        if tx_id in user_tx_colors: current_exon_color = user_tx_colors[tx_id]
        elif exon_override: current_exon_color = exon_override
        else: current_exon_color = tx_colors.get(tx_id, default_tx_color)

        # CDS Color
        base_color = tx_colors.get(tx_id, default_tx_color)
        if cds_override: current_cds_color = cds_override
        elif tx_id in user_tx_colors: current_cds_color = user_tx_colors[tx_id]
        elif base_color == default_tx_color: current_cds_color = default_cds_color
        else: current_cds_color = base_color

        for exon_start, exon_end in transcript_data['exons']:
            ax1.add_patch(patches.Rectangle((exon_start, i - 0.1), exon_end - exon_start, 0.2, color=current_exon_color))

        for cds_start, cds_end in transcript_data['cds']:
            ax1.add_patch(patches.Rectangle((cds_start, i - 0.2), cds_end - cds_start, 0.4, color=current_cds_color))

        # Features
        if tx_id not in model_features: continue
        Y_LOW, Y_UPP = i - 0.8, i - 0.22 
        coord_d = coordinate_adjustment(feat_names, Y_LOW, Y_UPP)

        for k_feat, feat in enumerate(feat_names):
            if feat not in model_features[tx_id]: continue
            y_low, y_upp = coord_d[feat]
            if not model_features[tx_id][feat]: continue
            
            plot_order = sorted(model_features[tx_id][feat], key=lambda k: sum(e - s for s, e in model_features[tx_id][feat][k]), reverse=True)
            plot_groups = overlaped_feature_indivisuals(model_features[tx_id][feat])
            plot_groups_cnt = [-1 for _ in range(len(plot_groups))]

            for ind in plot_order:
                idx = 0
                for g in range(len(plot_groups)):
                    if ind in plot_groups[g]:
                        plot_groups_cnt[g] += 1
                        idx = g
                        break
                y_low2 = y_low - 0.009 * plot_groups_cnt[idx]
                y_upp2 = y_upp - 0.009 * plot_groups_cnt[idx]
                crush_rate = 0.5
                full_height = y_upp2 - y_low2
                box_height = full_height * crush_rate
                y_center = (y_upp2 + y_low2) / 2
                y_start = y_center - (box_height / 2)

                for (start_k, end_k) in model_features[tx_id][feat][ind]:
                    ax1.add_patch(patches.Rectangle((start_k, y_start), end_k - start_k, box_height, color=feat_colors[k_feat]))
                leftmost_x = model_features[tx_id][feat][ind][0][0]
                rightmost_x = model_features[tx_id][feat][ind][-1][1]
                ax1.add_patch(patches.Rectangle((leftmost_x, y_start), rightmost_x - leftmost_x, box_height, linewidth=0.5, edgecolor=darken_color(feat_colors[k_feat]), facecolor="none"))

    if not hide_scale:
        data_width = MAX - MIN
        scale_size, scale_label = get_scale_unit(data_width)
        bar_y = -0.8
        text_y = bar_y - 0.05
        bar_end = MAX
        bar_start = MAX - scale_size
        ax1.plot([bar_start, bar_end], [bar_y, bar_y], color="#3B5894", linewidth=1, clip_on=False)
        cap_height = 0.1
        ax1.plot([bar_start, bar_start], [bar_y - cap_height/2, bar_y + cap_height/2], color="#3B5894", linewidth=1, clip_on=False)
        ax1.plot([bar_end, bar_end], [bar_y - cap_height/2, bar_y + cap_height/2], color="#3B5894", linewidth=1, clip_on=False)
        ax1.text((bar_start + bar_end) / 2, text_y, scale_label, ha='center', va='top', fontsize=8, color="#3B5894")
        if ci is not None:
            ax1.text(bar_end, bar_y + 0.1, "Note: Introns compressed", ha='right', va='bottom', fontsize=8, color="#3B5894", style='italic')

    space = int((MAX-MIN)/20)
    ax1.set_xlim(MIN - 2*space, MAX + space)
    ax1.set_ylim(-1.0, len(transcripts_data) - 0.5)

    if tss_mode:
        unique_tss = sorted(set(all_tss))
        ymin, ymax = -0.3, len(transcripts_data) - 0.8
        for x in unique_tss:
            ax1.vlines(x, ymin=ymin, ymax=ymax, colors=line_color, linestyles=':', linewidth=0.8, alpha=0.4)

    chrom_label = transcripts_data[0]["seq_region_name"]
    chrom_norm = chrom_label[3:] if chrom_label.lower().startswith("chr") else chrom_label
    ax1.set_xlabel(f'Chr{chrom_norm}', color=text_color)
    ax1.set_title(gene_name, color=text_color)
    ax1.set_yticks(y_positions)
    ax1.set_yticklabels([t['id'] for t in transcripts_data], color=text_color)

    N = len(transcripts_data) - 1
    if feat_names:
        coord = coordinate_adjustment(feat_names, Y_LOW=N - 0.8, Y_UPP=N - 0.22)
        for k_feat, feat in enumerate(feat_names):
            y_low, y_upp = coord[feat]
            ax1.text(MIN - 1.8 * space, (y_upp + y_low) / 2, feat, fontsize=9, color=darken_color(feat_colors[k_feat], factor=0.6), ha="left", va="center")

    for label in ax1.get_xticklabels():
        label.set_rotation(15)
        label.set_ha('right')
        label.set_color(text_color)
    return ax1


def get_diu_data(diu_file, transcripts_data, gene_symbol,
                 col_geneID="gene_id",
                 col_tx="tx_id",
                 col_prop_ctrl="prop_0",
                 col_prop_case="prop_1",
                 col_p="tx_p"):
    
    df = pd.read_csv(diu_file, sep="\t")
    
    # 必要なカラムの存在チェック (もしファイルにその名前がなければエラーになる)
    # ユーザー指定のフォーマットに準拠するため、ここは厳密に見る
    required_cols = [col_geneID, col_tx, col_prop_ctrl, col_prop_case, col_p]
    
    # 簡単なチェック (欠損があれば警告して戻る)
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"Error: The following required columns are missing in DIU file: {missing}")
        print(f"Expected: {required_cols}")
        print(f"Found: {list(df.columns)}")
        return transcripts_data

    # フィルタリング
    df_gene = df[df[col_geneID] == gene_symbol].copy()
    
    if df_gene.empty:
         print(f"Warning: No DIU data found for gene {gene_symbol}")
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
    if not order_list: return transcripts_data
    out = []
    for tx in order_list:
        for t in transcripts_data:
            if t['id'] == tx:
                out.append(t)
                break
    return out


def prepare_diu_ax2(transcripts_data, ax2, colors_meta, group0_label="Nontumor", group1_label="HCC"):
    tdata = transcripts_data[::-1]
    y_pos = np.arange(len(tdata))
    group0_vals = [t['group0_prop'] for t in tdata]
    group1_vals = [t['group1_prop'] for t in tdata]
    c_group0 = colors_meta['diu'].get('bar_group0', '#fdd9b5')
    c_group1 = colors_meta['diu'].get('bar_group1', '#d95f02')

    bar_h = 0.35
    ax2.barh(y_pos - bar_h/2, group1_vals, height=bar_h, label=group1_label, color=c_group1)
    ax2.barh(y_pos + bar_h/2, group0_vals, height=bar_h, label=group0_label, color=c_group0)

    max_usage = max(np.nanmax(group1_vals), np.nanmax(group0_vals))
    if not np.isfinite(max_usage) or max_usage <= 0: max_usage = 1.0

    ax2.set_xlim(0, max_usage * 1.05)
    ax2.set_ylim(-1.0, len(tdata) - 0.5)
    ax2.set_yticks([])
    ax2.tick_params(axis='y', length=0)
    from matplotlib.ticker import FuncFormatter
    ax2.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '0' if x == 0 else f'{x:.2f}'))
    ax2.set_xlabel("Mean usage")
    ax2.legend(fontsize=8, frameon=False)
    return ax2


def prepare_diu_ax3(transcripts_data, ax3, p_thresh=0.05, group0_label="0", group1_label="1"):
    tdata = transcripts_data[::-1]
    y_pos = np.arange(len(tdata))
    delta_vals = np.array([t['delta_prop'] for t in tdata], dtype=float)
    pvals      = np.array([t['pval_diu']   for t in tdata], dtype=float)
    
    valid_mask = (~np.isnan(pvals)) & (~np.isnan(delta_vals))
    sig_mask = (pvals < p_thresh) & valid_mask

    # --- log(p) 計算 ---
    with np.errstate(divide='ignore'):
        logp = -np.log10(pvals)
    
    max_logp_sig = 1.0
    if np.any(sig_mask):
        m = np.nanmax(logp[sig_mask])
        if np.isfinite(m) and m > 0: max_logp_sig = m

    # --- 色の決定 ---
    colors = []
    for du, pv, lp, sig in zip(delta_vals, pvals, logp, sig_mask):
        if (not sig) or (not np.isfinite(du)):
            colors.append((0.7, 0.7, 0.7)) # Gray
        else:
            intensity = lp / max_logp_sig
            intensity = max(0.0, min(1.0, intensity))
            if du >= 0: # Red
                colors.append((1.0, 0.8 - 0.5*intensity, 0.8 - 0.5*intensity))
            else:       # Blue
                colors.append((0.8 - 0.5*intensity, 0.8 - 0.5*intensity, 1.0))

    # --- バーの描画 ---
    ax3.barh(y_pos, delta_vals, height=0.6, color=colors, edgecolor="none")
    ax3.axvline(0, color="k", linewidth=0.8)

    # --- テキスト描画 (枠外配置) ---
    for y, d_val, p_val in zip(y_pos, delta_vals, pvals):
        # delta usage
        if np.isnan(d_val):
            d_text = "NA"
        else:
            d_text = f"{d_val:.2f}".replace('-', '–')

        # p value
        if np.isnan(p_val):
            raw_p_text = "NA"
        elif p_val < 0.001:
            raw_p_text = f"{p_val:.1e}".replace('-', '–')
        else:
            raw_p_text = f"{p_val:.3f}"

        if np.isnan(p_val):
            is_significant = False
        else:
            is_significant = p_val < p_thresh
        
        if is_significant:
            label_str = f"p={raw_p_text}\nΔ={d_text}"
        else:
            label_str = f"p={raw_p_text}\nΔ={d_text}"

        # 枠外に配置 (transform=ax3.get_yaxis_transform() 使用)
        ax3.text(1.05, y, label_str, 
                 ha='left', va='center', fontsize=10, color='black',
                 transform=ax3.get_yaxis_transform())

    # --- 軸範囲(xlim)の設定 ---
    max_abs = np.nanmax(np.abs(delta_vals[valid_mask])) if np.any(valid_mask) else 0.1
    if (not np.isfinite(max_abs)) or max_abs == 0: max_abs = 0.1
    
    pad = max_abs * 0.1
    ax3.set_xlim(-(max_abs + pad), max_abs + pad)

    # --- その他の整形 ---
    ax3.set_ylim(-1.0, len(tdata) - 0.5)
    
    from matplotlib.ticker import FuncFormatter
    ax3.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '0' if x == 0 else f'{x:.2f}'))
    
    ax3.set_yticks([])
    ax3.tick_params(axis='y', length=0)
    ax3.set_xlabel(f"Δusage ({group1_label} - {group0_label})")

    return ax3


# ---------------------------------------------------------
# Main Logic
# ---------------------------------------------------------
def isoespy_diu(gene_symbol, gtf_data, diu_tsv, meta_data, ci=None, group0_label="0", group1_label="1", p_thresh=0.05, tss_mode=False, output_file=None, dpi=300, hide_scale=False, show_features=False):
    config_meta, gtf_meta, colors_meta, query, feature_meta = parse_metadata(meta_data, gene_symbol)
    tx_colors = colors_meta['transcripts']
    transcripts_exon, transcripts_cds, tx_colors = get_isoform_model(gtf_data, gtf_meta, tx_colors, query, colors_meta)
    
    if show_features:
        model_features = get_feature_model(gtf_data, gtf_meta, query, feature_meta)
    else:
        model_features = {}

    transcripts_main = []
    transcripts_main = formatting_isoform_model(transcripts_main, transcripts_exon, annot="exon")
    transcripts_main = formatting_isoform_model(transcripts_main, transcripts_cds,  annot="cds")

    if not transcripts_main:
        print(f"No transcripts found for gene: {gene_symbol}")
        return

    x_min = min(transcripts_main, key=lambda x: x['start'])['start']
    x_max = max(transcripts_main, key=lambda x: x['end'])['end']

    ci_val = process_ci(ci)
    if ci_val is not None:
        copied_data = copy.deepcopy(transcripts_main)
        transcripts_main = moved_data_for_exons_cds(transcripts_main, ci_val)
        if show_features:
            model_features = moved_data_for_features(model_features, copied_data, ci_val)

    x_min_eff = min(transcripts_main, key=lambda x: x['start'])['start']
    x_max_eff = max(transcripts_main, key=lambda x: x['end'])['end']

    transcripts_main = get_diu_data(diu_tsv, transcripts_main, gene_symbol)
    transcripts_main = reorder(transcripts_main, config_meta.get('order'))
    
    if not transcripts_main:
        print("No transcripts left after filtering.")
        return

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(20,5), gridspec_kw={'width_ratios':[10,6,3]})

    ax1 = prepare_ax1_isoform_structure(transcripts_main, ax1, gene_symbol, tx_colors, colors_meta, tss_mode, model_features, feature_meta, ci_val, hide_scale)
    ax1 = prepare_ax1_xaxis(ax1, ci_val, x_min, x_max, x_min_eff, x_max_eff)
    ax2 = prepare_diu_ax2(transcripts_main, ax2, colors_meta, group0_label=group0_label, group1_label=group1_label)
    ax3 = prepare_diu_ax3(transcripts_main, ax3, p_thresh=p_thresh, group0_label=group0_label, group1_label=group1_label)

    plt.subplots_adjust(wspace=0.02)

    if output_file:
        try:
            fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
            print(f"Figure successfully saved to {output_file} with DPI={dpi}")
        except Exception as e:
            print(f"Error saving figure to file: {e}")
            plt.show()
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
    parser.add_argument("--hide_scale", action="store_true", help="Hide scale bar")
    parser.add_argument("-sf", "--show_features", action="store_true", help="Show genomic features (e.g. domains) in the plot")

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
        dpi=args.dpi,
        hide_scale=args.hide_scale,
        show_features=args.show_features
    )

if __name__ == "__main__":
    main()

