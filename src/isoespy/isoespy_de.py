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
from collections import defaultdict

try:
    from intronCompression import intronCompression
except ModuleNotFoundError:
    try:
        from isoespy.intronCompression import intronCompression
    except ModuleNotFoundError:
        def intronCompression(model, ci):
            return model

# helper for attributes
def get_attr_candidate(attr_dict, keys):
    """Try multiple keys to find a value in attr_dict."""
    for k in keys:
        if k in attr_dict:
            return attr_dict[k]
    return None


# metadata parser
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
        "de": {
            "box_group0": "#4daf4a",
            "box_group1": "#e41a1c"
        },
    }

    ctrl_trgt_table = dict()

    for line in lines:
        line = line.strip()
        if not line or line.startswith("!"):
            continue

        if line.startswith("[") and line.endswith("]"):
            current_section = line[1:-1]
            continue

        content = line.lstrip("#").strip() if line.startswith("#") else line

        if current_section == "sample":
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
                    if ":" in value:
                        val1, val2 = value.split(":", 1)
                        val1 = val1.strip()
                        val2 = [i.strip() for i in val2.split(",")]
                        for tx in val2:
                            if tx not in colors_meta["transcripts"]:
                                colors_meta["transcripts"][tx] = val1

        elif current_section == "gtf":
            if "=" in content:
                key, value = content.split("=", 1)
                gtf_meta[key.strip()] = value.strip()

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
                colors_meta["transcripts"][key.strip()] = value.strip()

        elif current_section == "colors.exon":
            if "=" in content:
                key, value = content.split("=", 1)
                if key.strip() == "color": colors_meta["exon"]["color"] = value.strip()
        
        elif current_section == "colors.cds":
            if "=" in content:
                key, value = content.split("=", 1)
                if key.strip() == "color": colors_meta["cds"]["color"] = value.strip()

        elif current_section == "colors.de":
            if "=" in content:
                key, value = content.split("=", 1)
                key = key.strip()
                if key in colors_meta["de"]:
                    colors_meta["de"][key] = value.strip()

    query["tx"][0] = gtf_meta.get("transcript_id")
    if query["tx"][1] == set():
        query["tx"][1] = None

    return sample_meta, config_meta, ctrl_trgt_table, gtf_meta, colors_meta, query, feature_meta


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


def get_isoform_model(gtf_file, gtf_meta, tx_colors, query, colors_meta=None):
    transcripts = {}
    transcripts_CDS = {}
    target_gene = query["gene"][1]
    target_tx = query["tx"][1]

    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            chrom = chrom.replace("chr","")
            
            attr_dict = {match.group(1): match.group(2) for match in re.finditer(r'(\S+)\s+"([^"]+)"', attributes)}
            
            transcript_id = attr_dict.get(gtf_meta['transcript_id'])
            line_gene = attr_dict.get(query["gene"][0])
            line_tx = attr_dict.get(query["tx"][0])
            
            if line_gene != target_gene: continue
            if target_tx is not None and line_tx not in target_tx: continue

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

    for tx in transcripts:
        transcripts[tx][0] = sorted(transcripts[tx][0], key=lambda x:x[0])
    for tx in transcripts_CDS:
        transcripts_CDS[tx][0] = sorted(transcripts_CDS[tx][0], key=lambda x:x[0])

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


# feature model helpers
def get_feature_model(gtf_file, gtf_meta, query, feature_meta):
    """Read positional functional features (GTF-based) with robust parsing."""
    print(f"DEBUG: Parsing features from {gtf_file}")
    ff_d = dict()

    target_gene = query["gene"][1]
    target_tx = query["tx"][1]
    
    # Candidate keys for gene name if primary key fails
    gene_keys = [query["gene"][0], 'gene_name', 'gene_symbol', 'gene_id', 'Name']
    gene_keys = [k for k in gene_keys if k]

    feature_count = 0
    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith("#"): continue
            fields = line.strip().split("\t")
            if len(fields) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            if feature in feature_meta:
                # Robust parsing 1: Regex
                attr_dict = {
                    match.group(1): match.group(2)
                    for match in re.finditer(r'(\S+)\s+"([^"]+)"', attributes)
                }
                # Robust parsing 2: Split by ; and =
                if not attr_dict:
                    for raw in attributes.strip().split(";"):
                        parts = raw.strip().split("=") 
                        if len(parts) == 2:
                            attr_dict[parts[0].strip()] = parts[1].strip().strip('"')

                transcript_id = get_attr_candidate(attr_dict, [gtf_meta.get("transcript_id"), "transcript_id", "transcript", "ID"])
                
                # Try finding gene name
                line_gene = get_attr_candidate(attr_dict, gene_keys)
                line_tx = get_attr_candidate(attr_dict, [query["tx"][0], "transcript_id", "transcript", "ID"])

                # Filtering logic
                match_gene = (line_gene == target_gene)
                match_tx = (target_tx is None) or (line_tx in target_tx)

                # Fallback: if transcript matches target list, accept even if gene name is missing
                if not match_gene and line_tx and target_tx and (line_tx in target_tx):
                    match_gene = True

                if not match_gene:
                    continue
                if not match_tx:
                    continue

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
                
                feature_count += 1

    for tx in ff_d:
        for feat in feature_meta:
            for ind in ff_d[tx][feat]:
                ff_d[tx][feat][ind] = sorted(ff_d[tx][feat][ind], key=lambda x: x[0])

    print(f"DEBUG: Total features accepted for gene {target_gene}: {feature_count}")
    return ff_d

def coordinate_adjustment(feat_names, Y_LOW, Y_UPP):
    out_d = dict()
    if feat_names == []:
        return out_d
    width = (Y_UPP - Y_LOW) / len(feat_names)
    for i in range(len(feat_names)):
        out_d[feat_names[i]] = (Y_UPP - width * (i + 1), Y_UPP - width * i)
    return out_d

def darken_color(color, factor=0.4):
    return tuple(max(0, c * factor) for c in color)

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
            transcripts_data[i]['det'] = [0, 1.0]

    return transcripts_data


def prepare_de_ax1(transcripts_data, ax1, gene_name, tx_colors, colors_meta, tss_mode, model_features, feature_meta, ci=None, hide_scale=False):
    # AX1: 各アイソフォームを視覚化
    y_positions = []
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
    default_tx_color = colors_meta.get("global", {}).get("default_tx", "#999999")
    default_cds_color = colors_meta.get("global", {}).get("default_cds", "#3b73b9")

    exon_override = colors_meta.get("exon", {}).get("color", None)
    cds_override = colors_meta.get("cds", {}).get("color", None)
    user_tx_colors = colors_meta.get("user_transcripts", {})

    for i, transcript_data in enumerate(transcripts_data):
        y_positions.append(i)
        start = transcript_data['start']
        end = transcript_data['end']
        tx_name = transcript_data['id']
        strand = transcript_data['strand']

        arrow_direction = "right" if strand == 1 else "left"
        y_pos = i
        
        ax1.annotate('', xy=(end, y_pos), xytext=(start, y_pos), arrowprops=dict(arrowstyle="-", color=line_color, lw=1))
        
        interval = max(1, (MAX-MIN)//50)
        x_positions = np.arange(start, end, interval)
        x_positions = x_positions[1:-1]
        
        marker_shape = ">" if arrow_direction == "right" else "<"
        for x in x_positions:
            ax1.scatter(x, y_pos, marker=marker_shape, color=line_color, s=10)

        if strand == 1:
            all_tss.append(start)
        else:
            all_tss.append(end)

        # Exon
        if tx_name in user_tx_colors:
            current_exon_color = user_tx_colors[tx_name]
        elif exon_override:
            current_exon_color = exon_override
        else:
            current_exon_color = tx_colors.get(tx_name, default_tx_color)

        # CDS
        base_color = tx_colors.get(tx_name, default_tx_color)
        if cds_override:
            current_cds_color = cds_override
        elif tx_name in user_tx_colors:
            current_cds_color = user_tx_colors[tx_name]
        elif base_color == default_tx_color:
            current_cds_color = default_cds_color
        else:
            current_cds_color = base_color

        for exon in transcript_data['exons']:
            exon_start = exon[0]
            exon_end = exon[1]
            ax1.add_patch(patches.Rectangle((exon_start, i - 0.1), exon_end - exon_start, 0.2, color=current_exon_color))

        for cds in transcript_data['cds']:
            cds_start = cds[0]
            cds_end = cds[1]
            ax1.add_patch(patches.Rectangle((cds_start, i - 0.2), cds_end - cds_start, 0.4, color=current_cds_color))

        if tx_name not in model_features:
            continue

        Y_LOW = i - 0.8
        Y_UPP = i - 0.22 
        coord_d = coordinate_adjustment(feat_names, Y_LOW, Y_UPP)

        for k_feat, feat in enumerate(feat_names):
            if feat not in model_features[tx_name]: continue
            y_low, y_upp = coord_d[feat]
            if not model_features[tx_name][feat]: continue
            
            plot_order = sorted(
                model_features[tx_name][feat],
                key=lambda k: sum(e - s for s, e in model_features[tx_name][feat][k]),
                reverse=True,
            )
            plot_groups = overlaped_feature_indivisuals(model_features[tx_name][feat])
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

                for (start_k, end_k) in model_features[tx_name][feat][ind]:
                    ax1.add_patch(
                        patches.Rectangle(
                            (start_k, y_start),
                            end_k - start_k,
                            box_height,
                            color=feat_colors[k_feat],
                        )
                    )

                leftmost_x = model_features[tx_name][feat][ind][0][0]
                rightmost_x = model_features[tx_name][feat][ind][-1][1]
                ax1.add_patch(
                    patches.Rectangle(
                        (leftmost_x, y_start),
                        rightmost_x - leftmost_x,
                        box_height,
                        linewidth=0.5,
                        edgecolor=darken_color(feat_colors[k_feat]),
                        facecolor="none",
                    )
                )

    if tss_mode:
        unique_tss = sorted(set(all_tss))
        ymin = -0.3
        ymax = len(transcripts_data) - 0.8
        for x in unique_tss:
            ax1.vlines(x, ymin=ymin, ymax=ymax, colors=line_color, linestyles=':', linewidth=0.8, alpha=0.4)

    # === スケールバー ===
    if not hide_scale:
        data_width = MAX - MIN
        scale_size, scale_label = get_scale_unit(data_width)
        bar_y = - 0.8
        text_y = bar_y - 0.05
        bar_end = MAX
        bar_start = MAX - scale_size
        ax1.plot([bar_start, bar_end], [bar_y, bar_y], color="#3B5894", linewidth=1, clip_on=False)
        cap_height = 0.1
        ax1.plot([bar_start, bar_start], [bar_y - cap_height/2, bar_y + cap_height/2], color="#3B5894", linewidth=1, clip_on=False)
        ax1.plot([bar_end, bar_end], [bar_y - cap_height/2, bar_y + cap_height/2], color="#3B5894", linewidth=1, clip_on=False)
        ax1.text((bar_start + bar_end) / 2, text_y, scale_label, ha='center', va='top', fontsize=8, color="#3B5894")
        if ci is not None:
            note_text = "Note: Introns compressed"
            ax1.text(bar_end, bar_y + 0.1, note_text, ha='right', va='bottom', fontsize=8, color="#3B5894", style='italic')
    # ========================

    space = int((MAX-MIN)/20)
    ax1.set_xlim(MIN - 2*space, MAX + space)
    ax1.set_ylim(-1.0, len(transcripts_data) - 0.5)
    
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

    N = len(transcripts_data) - 1
    if feat_names:
        coord = coordinate_adjustment(feat_names, Y_LOW=N - 0.8, Y_UPP=N - 0.22)
        for k_feat, feat in enumerate(feat_names):
            y_low, y_upp = coord[feat]
            ax1.text(
                MIN - 1.8 * space,
                (y_upp + y_low) / 2,
                feat,
                fontsize=9,
                color=darken_color(feat_colors[k_feat], factor=0.6),
                ha="left",
                va="center",
            )

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

    positions = []
    for i, name in enumerate(transcripts):
         positions.extend([i-0.12, i+0.12])

    boxplot = ax2.boxplot(data, positions=positions, vert=False, patch_artist=True, widths=0.2, showfliers=outliers)

    c_group0 = colors_meta['de'].get('box_group0', '#4daf4a')
    c_group1 = colors_meta['de'].get('box_group1', '#e41a1c')
    
    box_colors = [c_group1 if idx % 2 == 0 else c_group0 for idx in range(len(data))]
    
    for patch, color in zip(boxplot['boxes'], box_colors):
        patch.set_facecolor(color)

    for median in boxplot['medians']:
        median.set_color('black')

    ax2.set_ylim(-1.0, len(transcripts_data) - 0.5) 
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
    # AX3: ヒートマップ (q値) + q値/logFCのテキスト表示
    
    # === ax1と順序を合わせる ===
    transcripts_data = transcripts_data[::-1]
    
    # === 全体の描画範囲を等分する ===
    N = len(transcripts_data)
    Y_MIN = -1.0
    Y_MAX = N - 0.5
    TOTAL_HEIGHT = Y_MAX - Y_MIN
    UNIT_HEIGHT = TOTAL_HEIGHT / N 

    def color_mapping(nlogqval):
        if nlogqval > 20: return 1
        else: return 0.3 + nlogqval*0.035
    
    DET = list()
    for i in range(len(transcripts_data)):
        DET.append([transcripts_data[i]['id']] + transcripts_data[i]['det'])

    warm_cmap = plt.cm.Reds
    cold_cmap = plt.cm.Blues

    for i, val in enumerate(DET):
        logFC = float(DET[i][1]) # floatに変換
        qval = float(DET[i][2])  # floatに変換
        
        # 有意かどうか判定
        is_significant = qval < float(qval_threshold)

        # 色の決定
        if is_significant:
            if logFC > 0:
                color = warm_cmap(color_mapping(-np.log10(qval)))
            else:
                color = cold_cmap(color_mapping(-np.log10(qval)))
        else:
                color = (0.5, 0.5, 0.5, 0.5)

        # ヒートマップのバーを描画
        y_bottom = Y_MIN + (i * UNIT_HEIGHT)
        rect = patches.Rectangle(
            (0, y_bottom),
            0.1,  # バーの幅
            UNIT_HEIGHT,
            facecolor=color, 
            edgecolor='black', 
            linewidth=2
        )
        ax3.add_patch(rect)

        # === 【テキスト描画ロジック】 ===
        
        # 1. q値の文字列作成
        if qval < 0.001:
            # 非常に小さい値は指数表記 (例: q=1.2e-05)
            text_str = f"q={qval:.1e}\n"
        else:
            # 通常は小数第3位まで (例: q=0.052)
            text_str = f"q={qval:.3f}\n"
            
        # 2. 有意な場合のみ logFC を追加
        logfc_str = f"{logFC:.2f}".replace('-', '–')
        if is_significant:
            text_str += f"logFC={logfc_str}"
        else:
            text_str += f"logFC=Not Sig."

        # 3. テキストを描画
        ax3.text(
            0.15,  # バーの少し右
            y_bottom + (UNIT_HEIGHT / 2), 
            text_str,
            ha='left',      # 左揃え
            va='center',    # 上下中央揃え
            fontsize=10,     # 文字が長くなるので少し小さめに設定
            color='black'
        )
        # ==============================

    # カラーバーの描画設定（変更なし）
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
    ax3.set_ylim(Y_MIN, Y_MAX)
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
        model_features,
        feature_meta, # Added
        output_file=None,
        dpi=300,
        q_thresh=0.05,
        hide_scale=False):
    
    transcripts_data = reorder(transcripts_data, config_meta)
    
    if not transcripts_data:
        print("No transcripts found for plotting.")
        return

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(20, 5), gridspec_kw={'width_ratios': [10, 6, 2]})

    ax1 = prepare_de_ax1(transcripts_data, ax1, gene_name, tx_colors, colors_meta, tss_mode, model_features, feature_meta, ci, hide_scale)
    ax1 = prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff)

    ax2 = prepare_de_ax2(transcripts_data, ax2, config_meta, colors_meta, outliers, group0_label, group1_label)

    ax3 = prepare_de_ax3(transcripts_data, q_thresh, ax3)

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


def moved_data_for_exons_cds(main_data, ci):
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


def isoespy_de(gene, gtf_data, expression_data, det_data, meta_data, ci, outliers, group0_label, group1_label, tss_mode, output_file=None, dpi=300, q_thresh=0.05, hide_scale=False, show_features=False):
    sample_meta, config_meta, ctrl_trgt_table, gtf_meta, colors_meta, query, feature_meta = parse_metadata(meta_data, gene)

    tx_colors = colors_meta["transcripts"]
    model_exon, model_cds, tx_colors = get_isoform_model(gtf_data, gtf_meta, tx_colors, query, colors_meta)
    
    # Feature model (Conditional)
    if show_features:
        model_features = get_feature_model(gtf_data, gtf_meta, query, feature_meta)
    else:
        model_features = {}

    transcripts_main_data = list()
    transcripts_main_data = formatting_isoform_model(transcripts_main_data, model_exon, annot="exon")
    transcripts_main_data = formatting_isoform_model(transcripts_main_data, model_cds, annot="cds")
    
    if not transcripts_main_data:
        print(f"No transcripts found for gene: {gene}")
        return

    x_min = min(transcripts_main_data, key=lambda x: x['start'])['start']
    x_max = max(transcripts_main_data, key=lambda x: x['end'])['end']

    if ci != None:
        copied_data = copy.deepcopy(transcripts_main_data)
        transcripts_main_data = moved_data_for_exons_cds(transcripts_main_data, ci)
        if show_features:
            model_features = moved_data_for_features(model_features, copied_data, ci)
    
    x_min_eff = min(transcripts_main_data, key=lambda x: x['start'])['start']
    x_max_eff = max(transcripts_main_data, key=lambda x: x['end'])['end']

    transcripts_main_data = get_expression_data(expression_data, transcripts_main_data, ctrl_trgt_table, sample_meta)
    transcripts_main_data = get_det_data(det_data, transcripts_main_data)

    plot_isoespy_de(
        transcripts_main_data,
        config_meta,
        gene,
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
        model_features,
        feature_meta,
        output_file,
        dpi,
        q_thresh,
        hide_scale
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
    parser.add_argument("--hide_scale", action="store_true", help="Hide scale bar")
    parser.add_argument("-sf", "--show_features", action="store_true", help="Show genomic features (e.g. domains) in the plot")

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
    hide_scale = args.hide_scale
    show_features = args.show_features

    isoespy_de(gene, gtf_data, expression_data, det_data, meta_data, ci, outliers, group0_label, group1_label, tss_mode, output_file, dpi, q_thresh, hide_scale, show_features)


if __name__ == '__main__':
    main()

