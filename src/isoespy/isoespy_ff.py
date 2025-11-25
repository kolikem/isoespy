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
    try:
        from isoespy.intronCompression import intronCompression
    except ModuleNotFoundError:
        def intronCompression(model, ci): return model


# -----------------------------
# metadata parser (Standardized)
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
        "ff": {
            "continuous": {"cmap": "viridis"},
            "categorical": {"palette": "tab10", "colors": {}},
            "binary": {"0": "#cccccc", "1": "#000000"},
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
        # コメント記号 # を除去して key, value を取得
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
                    # legacy style: "#colors = #HEX: tx1,tx2" (support but hierarchical wins)
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
            # features usually defined like "#Protein: feature_id" or "Protein: feature_id"
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

        # ---------------- colors.ff.* ----------------
        elif current_section == "colors.ff.continuous":
            if "=" in content:
                key, value = content.split("=", 1)
                if key.strip() == "cmap": colors_meta["ff"]["continuous"]["cmap"] = value.strip()

        elif current_section == "colors.ff.categorical":
            if "=" in content:
                key, value = content.split("=", 1)
                key, value = key.strip(), value.strip()
                if key == "palette":
                    colors_meta["ff"]["categorical"]["palette"] = value
                else:
                    # Specific category colors
                    colors_meta["ff"]["categorical"]["colors"][key] = value

        elif current_section == "colors.ff.binary":
            if "=" in content:
                key, value = content.split("=", 1)
                key, value = key.strip(), value.strip()
                if key in ("0", "1"):
                    colors_meta["ff"]["binary"][key] = value

    # finalize query
    query["tx"][0] = gtf_meta.get("transcript_id")
    if query["tx"][1] == set():
        query["tx"][1] = None

    return sample_meta, config_meta, ctrl_trgt_table, gtf_meta, colors_meta, query, feature_meta


# -----------------------------
# isoform model
# -----------------------------
def get_isoform_model(gtf_file, gtf_meta, tx_colors, query, colors_meta=None):
    """Build exon/CDS models. Also complete transcript colors."""
    transcripts = {}
    transcripts_CDS = {}

    target_gene = query["gene"][1]
    target_tx = query["tx"][1]

    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            chrom = chrom.replace("chr", "")

            attr_dict = {
                match.group(1): match.group(2)
                for match in re.finditer(r'(\S+)\s+"([^"]+)"', attributes)
            }
            transcript_id = attr_dict.get(gtf_meta["transcript_id"])
            line_gene = attr_dict.get(query["gene"][0])
            line_tx = attr_dict.get(query["tx"][0])

            if line_gene != target_gene:
                continue
            if target_tx is not None and line_tx not in target_tx:
                continue

            if feature == gtf_meta["exon"]:
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = [[], 1 if strand == "+" else -1, chrom]
                transcripts[transcript_id][0].append((int(start), int(end)))

            if feature == gtf_meta.get("cds", "CDS"):
                if transcript_id not in transcripts_CDS:
                    transcripts_CDS[transcript_id] = [[], 1 if strand == "+" else -1, chrom]
                transcripts_CDS[transcript_id][0].append((int(start), int(end)))

    for tx in transcripts:
        transcripts[tx][0] = sorted(transcripts[tx][0], key=lambda x: x[0])
    for tx in transcripts_CDS:
        transcripts_CDS[tx][0] = sorted(transcripts_CDS[tx][0], key=lambda x: x[0])

    # complete colors
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


# -----------------------------
# feature model
# -----------------------------
def get_feature_model(gtf_file, gtf_meta, query, feature_meta):
    """Read positional functional features (GTF-based)."""
    ff_d = dict()

    target_gene = query["gene"][1]
    target_tx = query["tx"][1]

    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            if feature in feature_meta:
                attr_dict = {
                    match.group(1): match.group(2)
                    for match in re.finditer(r'(\S+)\s+"([^"]+)"', attributes)
                }
                transcript_id = attr_dict.get(gtf_meta["transcript_id"])

                line_gene = attr_dict.get(query["gene"][0])
                line_tx = attr_dict.get(query["tx"][0])
                if line_gene != target_gene:
                    continue
                if target_tx is not None and line_tx not in target_tx:
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

    for tx in ff_d:
        for feat in feature_meta:
            for ind in ff_d[tx][feat]:
                ff_d[tx][feat][ind] = sorted(ff_d[tx][feat][ind], key=lambda x: x[0])

    return ff_d


def formatting_isoform_model(transcripts_data, transcripts, annot):
    if annot == "exon":
        for transcript_id, exons in transcripts.items():
            if not exons[0]: continue
            start = min([i[0] for i in exons[0]])
            end = max(i[1] for i in exons[0])
            transcripts_data.append(
                {
                    "id": transcript_id,
                    "exons": exons[0],
                    "strand": exons[1],
                    "seq_region_name": exons[2],
                    "start": start,
                    "end": end,
                }
            )
    elif annot == "cds":
        for i in range(len(transcripts_data)):
            isomodel = transcripts_data[i]
            tx_id = isomodel["id"]
            if tx_id in transcripts:
                transcripts_data[i]["cds"] = transcripts[tx_id][0]
            else:
                transcripts_data[i]["cds"] = []
    return transcripts_data


# -----------------------------
# nonGTF annotations
# -----------------------------
def get_nonGTF_annotations(annotation_file, gene):
    if annotation_file is None:
        return {}

    nongtf_d = dict()
    with open(annotation_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or line == "":
                continue
            cols = line.split("\t")
            if len(cols) < 5: continue
            tx_id, gene_id, feat_name, feat_type, status = cols
            if gene != gene_id:
                continue
            if tx_id not in nongtf_d:
                nongtf_d[tx_id] = {}
            nongtf_d[tx_id][feat_name] = [feat_type, status]

    features_d = dict()
    for tx_id in nongtf_d:
        for feat in nongtf_d[tx_id]:
            feat_type = nongtf_d[tx_id][feat][0]
            if feat not in features_d:
                features_d[feat] = [feat_type, {}]
    for feat in features_d:
        features_d[feat][1] = {tx_id: None for tx_id in nongtf_d}
    for tx_id in nongtf_d:
        for feat in nongtf_d[tx_id]:
            status = nongtf_d[tx_id][feat][1]
            features_d[feat][1][tx_id] = status

    return features_d


def prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff):
    if ci is None:
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    else:
        ax1.set_xticks([x_min_eff, x_max_eff])
        ax1.set_xticklabels([str(x_min), str(x_max)])
    return ax1


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


# -----------------------------
# ax1 drawing
# -----------------------------
def prepare_ff_ax1(transcripts_data, model_features, ax1, gene_name, tx_colors, colors_meta, tss_mode):
    if not transcripts_data: return ax1
    
    MIN = min(transcripts_data, key=lambda x: x["start"])["start"]
    MAX = max(transcripts_data, key=lambda x: x["end"])["end"]
    # Reverse order for plotting (Top=First in list)
    transcripts_data = transcripts_data[::-1]

    y_positions = []
    
    if model_features == {}:
        feat_names = []
        feat_colors = []
    else:
        # Check first existing transcript in model_features to get keys
        first_key = next(iter(model_features), None)
        if first_key:
            feat_names = list(model_features[first_key].keys())
            feat_colors = sns.color_palette("husl", len(feat_names))
        else:
            feat_names = []
            feat_colors = []

    all_tss = []
    
    # global colors
    gcol = colors_meta.get("global", {})
    line_color = gcol.get("default_line", "gray")
    text_color = gcol.get("default_text", "black")
    default_tx_color = gcol.get("default_tx", "#B3C8CF")

    exon_override = colors_meta.get("exon", {}).get("color", None)
    cds_override = colors_meta.get("cds", {}).get("color", None)

    for i, transcript_data in enumerate(transcripts_data):
        y_positions.append(i)
        start = transcript_data["start"]
        end = transcript_data["end"]
        tx_name = transcript_data["id"]

        strand = transcript_data["strand"]
        arrow_direction = "right" if strand == 1 else "left"
        y_pos = i

        if strand == 1:
            all_tss.append(start)
        else:
            all_tss.append(end)

        interval = max(1, (MAX - MIN) // 50)
        x_positions = np.arange(start, end, interval)[1:-1]

        ax1.annotate(
            "",
            xy=(end, y_pos),
            xytext=(start, y_pos),
            arrowprops=dict(arrowstyle="-", color=line_color, lw=1),
        )

        if arrow_direction == "right":
            ax1.scatter(x_positions, [y_pos] * len(x_positions), marker=">", color=line_color, s=10)
        else:
            ax1.scatter(x_positions, [y_pos] * len(x_positions), marker="<", color=line_color, s=10)

        # exon
        current_exon_color = exon_override or tx_colors.get(tx_name, default_tx_color)
        for exon in transcript_data["exons"]:
            exon_start, exon_end = exon
            ax1.add_patch(
                patches.Rectangle(
                    (exon_start, i - 0.04),
                    exon_end - exon_start,
                    0.08,
                    color=current_exon_color,
                )
            )

        # CDS
        current_cds_color = cds_override or tx_colors.get(tx_name, default_tx_color)
        for cds in transcript_data["cds"]:
            cds_start, cds_end = cds
            ax1.add_patch(
                patches.Rectangle(
                    (cds_start, i - 0.08),
                    cds_end - cds_start,
                    0.16,
                    color=current_cds_color,
                )
            )

        # positional features
        if tx_name not in model_features:
            continue

        Y_LOW = i - 0.8
        Y_UPP = i - 0.11
        coord_d = coordinate_adjustment(feat_names, Y_LOW, Y_UPP)

        for k_feat, feat in enumerate(feat_names):
            if feat not in model_features[tx_name]: continue
            
            y_low, y_upp = coord_d[feat]
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
                crush_rate = 0.83

                for (start_k, end_k) in model_features[tx_name][feat][ind]:
                    ax1.add_patch(
                        patches.Rectangle(
                            (start_k, y_upp2 - (y_upp2 - y_low2) * crush_rate),
                            end_k - start_k,
                            (y_upp2 - y_low2) * crush_rate,
                            color=feat_colors[k_feat],
                        )
                    )

                leftmost_x = model_features[tx_name][feat][ind][0][0]
                rightmost_x = model_features[tx_name][feat][ind][-1][1]
                ax1.add_patch(
                    patches.Rectangle(
                        (leftmost_x, y_upp2 - (y_upp2 - y_low2) * crush_rate),
                        rightmost_x - leftmost_x,
                        (y_upp2 - y_low2) * crush_rate,
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
            ax1.vlines(
                x,
                ymin=ymin,
                ymax=ymax,
                colors=line_color,
                linestyles=":",
                linewidth=0.8,
                alpha=0.4,
            )

    space = int((MAX - MIN) / 20)
    ax1.set_xlim(MIN - 2 * space, MAX + space)
    ax1.set_ylim(-1.0, len(transcripts_data) - 0.5)
    
    chrom_label = transcripts_data[0]["seq_region_name"]
    chrom_norm = chrom_label.lower()
    if chrom_norm.startswith("chr"):
        chrom_norm = chrom_label[3:]
    else:
        chrom_norm = chrom_label

    ax1.set_xlabel(f'Chr{chrom_norm}', color=text_color)
    ax1.set_title(gene_name, color=text_color)
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    ax1.set_yticks(y_positions)
    transcripts = [i["id"] for i in transcripts_data]
    ax1.set_yticklabels(transcripts, color=text_color)

    # feature legend labels
    N = len(transcripts_data) - 1
    if feat_names:
        coord = coordinate_adjustment(feat_names, Y_LOW=N - 0.8, Y_UPP=N - 0.11)
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
        label.set_ha("right")
        label.set_color(text_color)

    return ax1


# -----------------------------
# axk drawing
# -----------------------------
def prepare_ff_axk(ax_l, transcripts_data, annotation_dataB, ff_colors_meta):
    # Same as before, logic for feature heatmap
    def is_float(v):
        try:
            float(v)
            return True
        except (ValueError, TypeError):
            return False

    def Color_binary(A):
        bmap = ff_colors_meta.get("binary", {})
        c0 = mcolors.to_rgba(bmap.get("0", "#cccccc"))
        c1 = mcolors.to_rgba(bmap.get("1", "#000000"))
        binary_colors = {"0": c0, "1": c1, None: "white"}
        return {key: binary_colors.get(str(value), "white") if value is not None else "white" for key, value in A.items()}

    def Color_categorical(A):
        unique_categories = sorted(set(v for v in A.values() if v is not None))

        cat_meta = ff_colors_meta.get("categorical", {})
        manual = cat_meta.get("colors", {})
        palette_name = cat_meta.get("palette", "tab10")

        category_colors = {cat: manual[cat] for cat in manual if cat in unique_categories}
        remaining = [cat for cat in unique_categories if cat not in category_colors]
        
        if remaining:
            auto_cols = None
            try:
                cmap = plt.get_cmap(palette_name, len(remaining))
                auto_cols = [cmap(i) for i in range(len(remaining))]
            except ValueError:
                auto_cols = None
            
            if auto_cols is None:
                try:
                    pal = sns.color_palette(palette_name, len(remaining))
                    auto_cols = [mcolors.to_rgba(c) for c in pal]
                except ValueError:
                    pal = sns.color_palette("tab10", len(remaining))
                    auto_cols = [mcolors.to_rgba(c) for c in pal]

            for cat, col in zip(remaining, auto_cols):
                category_colors[cat] = col

        category_colors[None] = "white"
        return {key: category_colors[value] for key, value in A.items()}

    def Color_continuous(A):
        valid_values = [float(v) for v in A.values() if v is not None and is_float(v)]
        if not valid_values:
            return {key: "white" for key in A}

        vmin, vmax = min(valid_values), max(valid_values)
        if vmin == vmax: vmax = vmin + 1

        cmap_name = ff_colors_meta.get("continuous", {}).get("cmap", "viridis")
        cmap = plt.get_cmap(cmap_name)

        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        return {
            key: cmap(norm(float(value)))
            if value is not None and is_float(value)
            else "white"
            for key, value in A.items()
        }

    # Transcripts are already reversed in plot_isoespy_ff before calling this, 
    # BUT prepare_ff_ax1 reverses them AGAIN for its own internal logic.
    # Here, transcripts_data is the ordered list (top to bottom of list).
    # We want top of list = top of plot (y=N-1).
    
    # Actually, plot_isoespy_ff calls: reorder -> (list is Order1, Order2...)
    # Then calls prepare_ff_ax1 -> reverses locally to plot.
    # Then calls prepare_ff_axk -> we need to plot Order1 at top.
    
    # In axk loop below: y = N - i - 0.9 - 1 ...
    # i=0 (Order1) -> y approx N. Correct.
    
    # fill missing txs
    for transcript_data in transcripts_data:
        tx_name = transcript_data["id"]
        for feat in annotation_dataB:
            if tx_name not in annotation_dataB[feat][1]:
                annotation_dataB[feat][1][tx_name] = None

    N = len(transcripts_data)
    features_l = list(annotation_dataB.keys())

    for i, transcript_data in enumerate(transcripts_data):
        tx_name = transcript_data["id"]

        for j, feat in enumerate(features_l):
            CATEGORY = annotation_dataB[feat][0]
            VALUES = annotation_dataB[feat][1]

            if CATEGORY == "binary":
                colormap = Color_binary(VALUES)
            elif CATEGORY == "categorical":
                colormap = Color_categorical(VALUES)
            elif CATEGORY == "continuous":
                colormap = Color_continuous(VALUES)
            else:
                colormap = {k: "white" for k in VALUES}

            hatch = "//" if VALUES[tx_name] is None else None
            rect = patches.Rectangle(
                (0.3, N - i - 0.9 - 1),
                0.4,
                1,
                facecolor=colormap.get(tx_name, "white"),
                edgecolor="black",
                hatch=hatch,
                linewidth=2,
            )
            ax_l[j].add_patch(rect)

            val_disp = VALUES[tx_name] if VALUES[tx_name] is not None else ""
            x = 0.3 + 0.4 / 2
            y = N - i - 0.9 - 1 + 1 / 2
            ax_l[j].text(
                x,
                y,
                val_disp,
                ha="center",
                va="center",
                rotation="vertical",
                fontsize=10,
            )

    for j, feat in enumerate(features_l):
        ax_l[j].set_title(feat)

    for ax in ax_l:
        ax.set_xlim(0, 1)
        ax.set_ylim(-1.0, N - 0.5)
        ax.axis("off")

    return ax_l


def reorder(transcripts_data, meta_data):
    if "order" not in meta_data:
        return transcripts_data
    tmp = []
    order = meta_data["order"]
    for tx in order:
        for i in transcripts_data:
            if i["id"] == tx:
                tmp.append(i)
                break
    return tmp


def plot_isoespy_ff(transcripts_data, model_features, config_meta, gene_name,
                    meta_data, tx_colors, colors_meta, ci,
                    x_min, x_max, x_min_eff, x_max_eff, annotation_dataB, tss_mode,
                    output_file=None, dpi=300):
    # 並び替え
    transcripts_data = reorder(transcripts_data, meta_data)
    
    if not transcripts_data:
        print("No transcripts found to plot.")
        return

    N = len(annotation_dataB.keys())
    if N == 0:
        fig, ax1 = plt.subplots(ncols=1, figsize=(20, 8))
        ax_l = []
    else:
        # width ratios: Main(20) : Heatmaps(1 each)
        fig, axes = plt.subplots(
            ncols=N + 1,
            figsize=(20, 8),
            gridspec_kw={"width_ratios": [20] + [1] * N},
        )
        if N == 0:
             ax1 = axes
             ax_l = []
        else:
             ax1 = axes[0]
             ax_l = axes[1:]

    ax1 = prepare_ff_ax1(transcripts_data, model_features, ax1, gene_name, tx_colors, colors_meta, tss_mode)
    ax1 = prepare_ax1_xaxis(ax1, ci, x_min, x_max, x_min_eff, x_max_eff)

    if N != 0:
        ax_l = prepare_ff_axk(ax_l, transcripts_data, annotation_dataB, colors_meta["ff"])

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
        tx_name = tx_data["id"]
        model[tx_name + "_exons"] = tx_data["exons"]
        model[tx_name + "_cds"] = tx_data["cds"]

    model_compressed = intronCompression(model, ci)

    startend_d = {}
    for tx_data in main_data:
        tx_name = tx_data["id"]
        if tx_name + "_exons" in model_compressed and model_compressed[tx_name + "_exons"]:
            start = model_compressed[tx_name + "_exons"][0][0]
            end = model_compressed[tx_name + "_exons"][-1][1]
            startend_d[tx_name] = {"start": start, "end": end}
        else:
             startend_d[tx_name] = {"start": tx_data['start'], "end": tx_data['end']}

    for tx_data in main_data:
        tx_name = tx_data["id"]
        tx_data["exons"] = model_compressed[tx_name + "_exons"]
        tx_data["cds"] = model_compressed[tx_name + "_cds"]
        tx_data["start"] = startend_d[tx_name]["start"]
        tx_data["end"] = startend_d[tx_name]["end"]

    return main_data


def moved_data_for_features(model_features, main_data, ci):
    model = {}
    # ベースのExon/CDSを登録
    for tx_data in main_data:
        tx_name = tx_data["id"]
        model[tx_name + "_exons"] = tx_data["exons"]
        model[tx_name + "_cds"] = tx_data["cds"]

    # 特徴量の座標を登録
    for tx_name in model_features:
        for feat_name in model_features[tx_name]:
            for id_name in model_features[tx_name][feat_name]:
                model[tx_name + "_" + feat_name + "_" + id_name] = model_features[tx_name][feat_name][id_name]

    model_compressed = intronCompression(model, ci)

    # 変換後の座標を戻す
    for tx_name in model_features:
        for feat_name in model_features[tx_name]:
            for id_name in model_features[tx_name][feat_name]:
                key = tx_name + "_" + feat_name + "_" + id_name
                if key in model_compressed:
                    model_features[tx_name][feat_name][id_name] = model_compressed[key]

    return model_features


def isoespy_ff(gene, gtf_data, meta_data, ci, annotation_file, tss_mode, output_file=None, dpi=300):
    sample_meta, config_meta, ctrl_trgt_table, gtf_meta, colors_meta, query, feature_meta = parse_metadata(meta_data, gene)

    # isoform model + colors completion
    tx_colors = colors_meta["transcripts"]
    model_exon, model_cds, tx_colors = get_isoform_model(gtf_data, gtf_meta, tx_colors, query, colors_meta)

    # feature model (Genomic Feature)
    model_features = get_feature_model(gtf_data, gtf_meta, query, feature_meta)

    # formatted isoform model
    transcripts_main_data = []
    transcripts_main_data = formatting_isoform_model(transcripts_main_data, model_exon, annot="exon")
    transcripts_main_data = formatting_isoform_model(transcripts_main_data, model_cds, annot="cds")
    
    if not transcripts_main_data:
        print(f"No transcripts found for gene: {gene}")
        return

    x_min = min(transcripts_main_data, key=lambda x: x["start"])["start"]
    x_max = max(transcripts_main_data, key=lambda x: x["end"])["end"]

    if ci is not None:
        copied_data = copy.deepcopy(transcripts_main_data)
        transcripts_main_data = moved_data_for_exons_cds(transcripts_main_data, ci)
        # copyした元データを使って feature も圧縮対応
        model_features = moved_data_for_features(model_features, copied_data, ci)

    x_min_eff = min(transcripts_main_data, key=lambda x: x["start"])["start"]
    x_max_eff = max(transcripts_main_data, key=lambda x: x["end"])["end"]

    # Functional annotation data (heatmap)
    annotation_dataB = get_nonGTF_annotations(annotation_file, gene)

    plot_isoespy_ff(
        transcripts_main_data,
        model_features,
        config_meta,
        gene,
        config_meta,
        tx_colors,
        colors_meta,
        ci,
        x_min,
        x_max,
        x_min_eff,
        x_max_eff,
        annotation_dataB,
        tss_mode,
        output_file,
        dpi
    )


def process_ci(ci):
    if ci is None:
        return None
    if isinstance(ci, str):
        try:
            ci = float(ci)
        except ValueError:
            warnings.warn(f"Warning: ci should be a float/int, but received string '{ci}' that cannot be converted.")
            return None
    else:
        try:
            ci = float(ci)
        except (ValueError, TypeError):
            warnings.warn(f"Warning: ci should be a float/int, but received '{ci}' that cannot be converted.")
            return None
    return ci


def main(args=None):
    parser = argparse.ArgumentParser(description="isoespy_ff()")
    parser.add_argument("-gene", "--gene_name", required=True, type=str, default=None, help="Gene name")
    parser.add_argument("-gtf", "--gtf_data", required=True, type=str, default=None, help="GTF file")
    parser.add_argument("-meta", "--meta_data", required=True, type=str, default=None, help="metadata")
    parser.add_argument("-ci", "--compress_introns", default=None, help="intron compression parameter")
    parser.add_argument("-a", "--annotation", type=str, default=None, help="non-gtf transcript annotation file")
    parser.add_argument("-tss", "--tss_line", action="store_true", help="Show TSS support lines")
    parser.add_argument("-o", "--output_file", type=str, default=None, help="Save figure to file (e.g., plot.pdf or plot.png)")
    parser.add_argument("--dpi", type=int, default=300, help="Resolution for image formats (e.g., PNG) in DPI")


    args = parser.parse_args()
    gene = args.gene_name
    gtf_data = args.gtf_data
    meta_data = args.meta_data
    ci = process_ci(args.compress_introns)
    annotation_data = args.annotation
    tss_mode = args.tss_line
    output_file = args.output_file
    dpi = args.dpi

    isoespy_ff(gene, gtf_data, meta_data, ci, annotation_data, tss_mode, output_file, dpi)


if __name__ == "__main__":
    main()
