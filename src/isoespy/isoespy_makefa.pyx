""" GTF --> fasta
    python3 isoespy_makefa.py --gtf <GTF file> --genome <genome ref> --meta <metadata> --feature <exon/CDS> --type <nucleotide/amino_acid>
"""
import argparse
import pysam

def get_target_seq(genome_ref, chrom, start, end):
   # 参照ゲノムファイルを開く（.fai インデックスが必要）
   fasta_path = genome_ref
   fasta = pysam.FastaFile(fasta_path)
   sequence = fasta.fetch(chrom, start-1, end)  # fetchは0-based
   fasta.close()
   return sequence


def parse_meta(meta_file):
    transcript_id = "transcript_id"
    cds = "CDS"
    exon = "exon"

    gtf_field = 0
    with open(meta_file) as f:
        for line in f:
            if line.startswith("[gtf]"):
                gtf_field = 1
                continue
            elif line.startswith("["):
                gtf_field = 0
            if gtf_field == 1:
                if line.strip() == "":
                    continue
                feature0 = line.split("=")[0].strip()[1:]
                feature1 = line.split("=")[1].strip()
                if feature0 == "transcript_id":
                    transcript_id = feature1
                if feature0 == "cds":
                    cds = feature1
                if feature0 == "exon":
                    exon = feature1

    return {"transcript_id": transcript_id, "exon": exon, "cds": cds}


def parse_gtf(gtf_file, meta_data, target_feature):
    """ GTFファイルをパースし転写産物ごとのfeature情報を取得 """
    gtf_data = {}
    
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):  # コメント行はスキップ
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, frame, attributes = cols
            start, end = int(start), int(end)
            
            # transcript_idを取得
            attr_dict = {k.strip(): v.strip('"') for k, v in 
                         (item.split() for item in attributes.split(';') if item)}
            transcript_id = attr_dict.get(meta_data["transcript_id"])

            if transcript_id:
                if feature == target_feature:
                    if transcript_id not in gtf_data:
                        gtf_data[transcript_id] = {"chrom": chrom, "strand": strand, target_feature: []}
                    gtf_data[transcript_id][target_feature].append((start, end))


    # 座標順への並び替え
    for tx in gtf_data:
        gtf_data[tx][target_feature] = sorted(gtf_data[tx][target_feature], key=lambda x: x[0])

    return gtf_data


def reverse_complement(seq):
    seq = seq.upper()
    complement = str.maketrans("ATGC", "TACG")
    return seq.translate(complement)[::-1]


def NUCSEQ_for_all_tx(gtf_data, genome_ref, feature):
    NUCSEQs = {}
    for tx_i in gtf_data:
        seq = ""
        chrom = gtf_data[tx_i]["chrom"]
        for cds_j in gtf_data[tx_i][feature]:
            seq += get_target_seq(genome_ref, chrom, cds_j[0], cds_j[1])
        NUCSEQs[tx_i] = seq.upper()

    # - ｽﾄﾗﾝﾄﾞはﾘﾊﾞｰｽｺﾝﾌﾟﾘﾒﾝﾄ
    for tx_i in gtf_data:
        if gtf_data[tx_i]["strand"] == "-":
            NUCSEQs[tx_i] = reverse_complement(NUCSEQs[tx_i])
    return NUCSEQs


def codon_table():
    # コドン表
    codon_table = {
    "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
    "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
    "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
    "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
    "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
    "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
    "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
    "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
    "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
    "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
    "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
    "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
    "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
    "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
    "TAC":"Y", "TAT":"Y", "TAA":"*", "TAG":"*",
    "TGC":"C", "TGT":"C", "TGA":"*", "TGG":"W"
    }
    return codon_table


def get_longest_orf_index(sequences):
    """
    アミノ酸配列のリストから、最初のストップコドン（*）までの
    読み枠が最も長い配列のインデックスを返す関数
    """
    orf_lengths = [seq.split('*')[0] for seq in sequences]
    max_index = max(range(len(orf_lengths)), key=lambda i: len(orf_lengths[i]))
    return max_index


def get_longest_orf_sequence(sequences):
    """
    アミノ酸配列のリストから、最初のストップコドン（*）までの
    読み枠が最も長い配列を返す関数
    """
    orfs = [seq.split('*')[0] for seq in sequences]
    longest_orf = max(orfs, key=len)
    return longest_orf


def AASEQ_for_all_tx(gtf_data, NUCSEQs, feature):
    if feature.upper() == "CDS":
        AASEQs = {}
        Frames = {}
        for tx_i in gtf_data:
            aa_seq_l = [[], [], []]    # frame 0, 1, 2
            # frame = 0
            nucseq_i = NUCSEQs[tx_i]
            for j in range(0, len(nucseq_i)-2, 3):
                codon = nucseq_i[j:j+3]
                aa_seq_l[0].append(codon_table().get(codon, "X")) # コドン表にない --> X
            # frame = 1
            nucseq_i = NUCSEQs[tx_i][1:]
            for j in range(0, len(nucseq_i)-2, 3):
                codon = nucseq_i[j:j+3]
                aa_seq_l[1].append(codon_table().get(codon, "X"))
            # frame = 2
            nucseq_i = NUCSEQs[tx_i][2:]
            for j in range(0, len(nucseq_i)-2, 3):
                codon = nucseq_i[j:j+3]
                aa_seq_l[2].append(codon_table().get(codon, "X"))

            # 3 seq
            aa_seq_l = ["".join(aa_seq_l[0]), "".join(aa_seq_l[1]), "".join(aa_seq_l[2])]

            # determine frame
            if not "*" in aa_seq_l[0]:
                AASEQs[tx_i] = aa_seq_l[0]
                Frames[tx_i] = "0"
            elif not "*" in aa_seq_l[1]:
                AASEQs[tx_i] = aa_seq_l[1]
                Frames[tx_i] = "1"
            elif not "*" in aa_seq_l[2]:
                AASEQs[tx_i] = aa_seq_l[2]
                Frames[tx_i] = "2"
            else:
                AASEQs[tx_i] = get_longest_orf_sequence(aa_seq_l)
                longest_idx = get_longest_orf_index(aa_seq_l)
                Frames[tx_i] = str(longest_idx)
                print("Caution!!   ", tx_i, " cannot be translated without PTC. Adopted the largest frame = ", longest_idx)

    else:
        AASEQs = {}
        Frames = None
        for tx_i in gtf_data:
            nucseq_i = NUCSEQs[tx_i]
            aa_seq = []
            for j in range(0, len(nucseq_i)-2, 3):
                codon = nucseq_i[j:j+3]
                aa_seq.append(codon_table().get(codon, "X"))   # コドン表にないコドンはXとする
            AASEQs[tx_i] = "".join(aa_seq)
    return AASEQs, Frames


def output(SEQs, outfile):
    with open(outfile, mode="w") as f:
        for tx in SEQs:
            f.write(">"+tx+"\n")
            f.write(SEQs[tx]+"\n")


def output_frames(frames, outfile):
    with open(outfile, mode="w") as f:
        for tx in frames:
            f.write(tx + "\t" + frames[tx] + "\n")


def main(args=None):
    # parse
    parser = argparse.ArgumentParser(description='isoespy_makefa()')
    parser.add_argument('--gtf', required=True, help='GTF file')
    parser.add_argument('--genome', required=True, help='genome reference')
    parser.add_argument('--meta', required=True, help='metadata file')
    parser.add_argument('--feature', required=True, help='A feature that needs to be transformed to fasta (exon / CDS)')
    parser.add_argument('--type', required=True, choices=["nucleotide", "amino_acid"], help='nucleotide or amino_acid')

    args = parser.parse_args()
    gtf_file = args.gtf
    genome_ref = args.genome
    meta_file = args.meta
    feature = args.feature
    seqtype = args.type

    meta_data = parse_meta(meta_file)
    gtf_data = parse_gtf(gtf_file, meta_data, feature)
    NUCSEQs = NUCSEQ_for_all_tx(gtf_data, genome_ref, feature)
    if seqtype == "nucleotide":
        # 塩基配列を出力
        output(NUCSEQs, gtf_file+f".{feature}.nt.fa")
    else:
        # アミノ酸配列を出力
        AASEQs, Frames = AASEQ_for_all_tx(gtf_data, NUCSEQs, feature)
        output(AASEQs, gtf_file+f".{feature}.aa.fa")
        if Frames != None:
            # AA CDS の場合のみ frame ファイルを出力
            output_frames(Frames, gtf_file+f".{feature}.aa.frame")

if __name__ == "__main__":
    main()

