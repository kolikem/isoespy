import argparse
import csv
import sys
import re
from collections import OrderedDict

def parse_args():
    p = argparse.ArgumentParser(
        description="Create DRIMSeq input table (gene_id, feature_id, counts...) from TSV-GTF and counts matrix."
    )
    p.add_argument("--gtf-tsv", required=True, help="GTF as TSV (attributes in 9th column)")
    p.add_argument("--counts", required=True, help="Counts matrix TSV (col1=tx, col2+=samples)")
    p.add_argument("--gene-key", default="gene_symbol", help="Attribute key for gene id (default: gene_symbol)")
    p.add_argument("--tx-key", default="transcript_id", help="Attribute key for transcript id (default: transcript_id)")
    p.add_argument("--output", default="drimseq_counts.tsv", help="Output TSV (default: drimseq_counts.tsv)")
    return p.parse_args()

def parse_attributes(attr_field):
    """
      - key "value"; key "value";
      - key=value;key=value;
      - key value; のような形式
    """
    attr = {}
    if attr_field is None:
        return attr
    s = attr_field.strip().strip(";")
    if not s:
        return attr

    # ; で分割（連続セミコロンや末尾空要素を無視）
    parts = re.split(r';\s*', s)
    for part in parts:
        part = part.strip()
        if not part:
            continue
        key = None
        val = None

        if '=' in part:
            key, val = part.split('=', 1)
        else:
            # key "value" あるいは key value を拾う
            m = re.match(r'(\S+)\s+"?(.*?)"?$', part)
            if m:
                key, val = m.group(1), m.group(2)

        if key is not None:
            key = key.strip()
            if val is not None:
                val = val.strip().strip('"')
            attr[key] = val
    return attr

def build_tx_to_gene_map(gtf_tsv_path, gene_key, tx_key):
    """
    第9列 attributes から gene_key と tx_key を取り出し、
    転写産物ID -> 遺伝子ID の辞書を作る。
    同一転写産物IDが複数行に出ても最初の対応を採用（2回目以降は無視）。
    """
    tx2gene = OrderedDict()
    line_no = 0
    with open(gtf_tsv_path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            line_no += 1
            if not row or len(row) < 9:
                # 不完全行はスキップ
                continue
            attributes_field = row[8]
            attrs = parse_attributes(attributes_field)

            tx = attrs.get(tx_key)
            gene = attrs.get(gene_key)

            if tx is None or gene is None:
                # 必要キーが無い行はスキップ（構造上、exon行など全てに両方あるとは限らないため）
                continue

            if tx not in tx2gene:
                tx2gene[tx] = gene
            else:
                # 2回目以降は飛ばす（ただし gene が不一致なら警告）
                if tx2gene[tx] != gene:
                    print(
                        f"[WARN] Transcript {tx!r} already mapped to {tx2gene[tx]!r}, "
                        f"but found another gene {gene!r} at line {line_no}. Keeping the first mapping.",
                        file=sys.stderr
                    )
    if not tx2gene:
        print("[ERROR] No (transcript_id, gene_id) pairs were extracted from attributes. "
              f"Check that 9th column contains attributes and keys '{tx_key}', '{gene_key}'.",
              file=sys.stderr)
        sys.exit(1)
    return tx2gene

def read_counts(counts_path):
    """
    counts TSV を読み込み：
      - 1行目: ヘッダー（col1=tx、col2+=サンプル名）
      - 2行目以降: tx, sample1_count, sample2_count, ...
    戻り値:
      header_samples: [sample1, sample2, ...]
      rows: [(tx, [c1, c2, ...]), ...]  ※入力の順序維持
    """
    rows = []
    with open(counts_path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            print("[ERROR] Counts file is empty.", file=sys.stderr)
            sys.exit(1)
        if len(header) < 2:
            print("[ERROR] Counts header must have at least two columns: tx and >=1 sample.", file=sys.stderr)
            sys.exit(1)

        # 1列目は tx、2列目以降はサンプル名
        sample_headers = header[1:]

        for i, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) < 2:
                print(f"[ERROR] Counts row {i} has fewer than 2 columns.", file=sys.stderr)
                sys.exit(1)
            tx = row[0]
            counts = row[1:]
            rows.append((tx, counts))

    if not rows:
        print("[ERROR] Counts file has no data rows.", file=sys.stderr)
        sys.exit(1)

    return sample_headers, rows

def ensure_sets_match(tx2gene, counts_rows):
    set_A = set(tx2gene.keys())
    set_B = set(tx for tx, _ in counts_rows)

    if set_A != set_B:
        only_in_A = sorted(set_A - set_B)
        only_in_B = sorted(set_B - set_A)
        msg_lines = ["[ERROR] Transcript ID sets do not match between GTF(attributes) and counts."]
        if only_in_A:
            msg_lines.append(f"  - Present only in GTF (A): {len(only_in_A)} examples: {only_in_A[:10]}")
        if only_in_B:
            msg_lines.append(f"  - Present only in counts (B): {len(only_in_B)} examples: {only_in_B[:10]}")
        msg_lines.append("This indicates some transcripts appear in only A or only B. "
                         "Please resolve the mismatch in your inputs.")
        print("\n".join(msg_lines), file=sys.stderr)
        sys.exit(1)

def write_output(output_path, sample_headers, counts_rows, tx2gene):
    """
    DRIMSeq 形式のTSV:
      gene_id  feature_id  sample1  sample2  ...
      (fileBの行順を維持)
    """
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        header = ["gene_id", "feature_id"] + sample_headers
        writer.writerow(header)

        for tx, counts in counts_rows:
            gene = tx2gene.get(tx)
            if gene is None:
                # 通常ここには来ない（事前一致チェック済み）
                print(f"[WARN] tx {tx!r} not found in mapping; skipping.", file=sys.stderr)
                continue
            writer.writerow([gene, tx] + counts)

def main():
    args = parse_args()

    tx2gene = build_tx_to_gene_map(args.gtf_tsv, args.gene_key, args.tx_key)
    sample_headers, counts_rows = read_counts(args.counts)
    ensure_sets_match(tx2gene, counts_rows)
    write_output(args.output, sample_headers, counts_rows, tx2gene)

    print(f"[OK] Wrote DRIMSeq counts table: {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()

