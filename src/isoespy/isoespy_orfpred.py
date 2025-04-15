""" CPATとTransDecoderを実行し、unannotated転写産物のCDS/ORF予測をGTFファイルに組み入れるスクリプト
    CPATのnon-coding予測転写産物以外について、BEDファイル(TransDecoder)のCDS開始・終了情報をGTFのexon情報をもとにゲノム座標に変換し、CDS情報がGTFに存在しない場合は追加する処理
"""

import argparse
import os
import subprocess
from datetime import datetime
try:
    from isoespy_makegtf import isoespy_makegtf
except ModuleNotFoundError:
    from isoespy.isoespy_makegtf import isoespy_makegtf

# BED
#				  total_length																start	end
# ADGRF4_novel_tx_1       0       3073    	ID=ADGRF4_novel_tx_1.p1;GENE.ADGRF4_novel_tx_1~~ADGRF4_novel_tx_1.p1;ORF_type:complete_(+),score=87.78  0       +       620     2264    0       1       3073    0
# 1-based座標系では、621塩基めから2264塩基めまでがCDS予測.

#	## transcript-based annotation format
#       #Transcript             from    to      coordinate      feature                 source          attributes
#       ADGRF4_novel_tx_1       252     500     amino_acid      protein_domain          pfam            domain_id "PF00002.29"; domain_name "7tm_2"; evalue "3.2e-33";
#       ADGRF4_novel_tx_1       199     242     amino_acid      protein_domain          pfam            domain_id "PF01825.27"; domain_name "GPS"; evalue; "4e-10";
#       ENST00000283303.3       252     500     amino_acid      protein_domain          pfam            domain_id "PF00002.29"; domain_name "7tm_2"; evalue "5.5e-33";
#       ENST00000283303.3       347     390     amino_acid      protein_domain          pfam            domain_id "PF01825.27"; domain_name "GPS"; evalue "5.4e-10";


def run_cpat(workdir, hexamer, model, fasta, outprefix):
    outprefix = os.path.basename(outprefix)

    # 作業ディレクトリが存在しない場合は作成
    os.makedirs(workdir, exist_ok=True)

    # エラー出力ファイルのパス
    now = datetime.now()
    timestamp = now.strftime("%y%m%d_%H%M%S")
    stderr_path = os.path.join(workdir, f"cpat_stderr_{timestamp}.txt")

    # CPATコマンド
    cmd = [
        "cpat",
        "-x", hexamer,
        "-d", model,
        "-g", fasta,
        "-o", outprefix
    ]

    # subprocessで実行、stderrをファイルにリダイレクト
    with open(stderr_path, "w") as err_log:
        subprocess.run(cmd, cwd=workdir, stderr=err_log)


def run_transdecoder(transdecoder_path, fasta, workdir):
    # ==== longorf ====
    exe_longorf = "TransDecoder.LongOrfs"
    if transdecoder_path:
        if os.path.basename(transdecoder_path) == exe_longorf:
            # 実行ファイルまでのフルパス指定のケース
            longorfs_cmd = transdecoder_path
        else:
            # ディレクトリ指定のケース
            longorfs_cmd = os.path.join(transdecoder_path, exe_longorf)
    else:
        longorfs_cmd = exe_longorf

    # LongOrfsコマンド
    os.makedirs(workdir, exist_ok=True)
    cmd = [
           longorfs_cmd,
           "-S",
           "-t", fasta
    ]

    now = datetime.now()
    timestamp = now.strftime("%y%m%d_%H%M%S")
    stderr_path = os.path.join(workdir, f"transdecoder_longorfs_stderr_{timestamp}.txt")

    # subprocessで実行、stderrをファイルにリダイレクト
    with open(stderr_path, "w") as err_log:
        subprocess.run(cmd, cwd=workdir, stderr=err_log)

    # ==== predict ====
    exe_predict = "TransDecoder.Predict"
    if transdecoder_path:
        if os.path.basename(transdecoder_path) == exe_predict:
            predict_cmd = transdecoder_path
        else:
            predict_cmd = os.path.join(transdecoder_path, exe_predict)
    else:
        predict_cmd = exe_predict

    # Predictコマンド
    cmd = [
           predict_cmd, 
           "--single_best_only",
           "-t", fasta
    ]

    now = datetime.now()
    timestamp = now.strftime("%y%m%d_%H%M%S")
    stderr_path = os.path.join(workdir, f"transdecoder_predict_stderr_{timestamp}.txt")

    # subprocessで実行、stderrをファイルにリダイレクト
    with open(stderr_path, "w") as err_log:
        subprocess.run(cmd, cwd=workdir, stderr=err_log)


def parse_bed(bed_file):
    """BEDファイルをパースし、転写産物ごとのCDS情報を取得"""
    bed_info = {}
    with open(bed_file) as f:
        for line in f:
            if line.startswith("track name="):
                continue
            cols = line.strip().split('\t')
            transcript_id = cols[0]
            cds_start = int(cols[6]) + 1   # 転写産物上のCDS開始    *1-based-inclusive に直す.
            cds_end = int(cols[7])         # 転写産物上のCDS終了
            strand = cols[5]               # ストランド (転写産物上におけるCDSの向き)
            total_len = int(cols[2])       # 転写産物全長

            # 終始コドンの削除 (GTFのCDS featureに含めないため)   *TransdecoderのCDS位置には終始コドンが含まれる
            cds_end -= 3
            
            bed_info[transcript_id] = {"cds_start": cds_start, "cds_end": cds_end, "total_len": total_len, "strand": strand}
    return bed_info


def parse_cpat(cpat_file):
    noncoding_s = set()
    with open(cpat_file) as f:
        for line in f:
            line = line.strip()
            if not line == "":
                noncoding_s.add(line)
    return noncoding_s


def parse_meta(meta_file):
    gtf_field = 0
    transcript_id = "transcript_id"
    cds = "CDS"
    with open(meta_file) as f:
        for line in f:
            if line.startswith("[gtf]"):
                gtf_field = 1
                continue
            elif line.startswith("["):
                gtf_field = 0
                continue
            if gtf_field == 1:
                if line.startswith("#transcript_id"):
                    transcript_id = line.strip().split("=")[1].strip()
                if line.startswith("#cds"):
                    cds = line.strip().split("=")[1].strip()
    return {"transcript_id": transcript_id, "cds": cds}


def parse_gtf(gtf_file, meta_d):
    cds_annotated_tx = dict()    # CDSがannotationされている転写産物は1, unannotatedなら0
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
            transcript_id = attr_dict.get(meta_d["transcript_id"])

            # cds_annotated_tx を更新
            if not transcript_id in cds_annotated_tx:
                cds_annotated_tx[transcript_id] = 0   # 0: CDS行ﾅｼ、1: CDS行ｱﾘ
            if feature == meta_d["cds"]:
                cds_annotated_tx[transcript_id] = 1
    return cds_annotated_tx


def add_cds_to_gtf_from_cpat_transdecoder(workdir, trans_outfile, cpat_outfile, gtf_file, meta_file):
    # CPATとTransDecoderの結果 --> isoespy_makegtf.py 入力形式に変換
    transdecoder_result = parse_bed(trans_outfile)    # TransDecoder の結果辞書
    cpat_result = parse_cpat(cpat_outfile)            # non-coding 判定の転写産物名セット

    meta_d = parse_meta(meta_file)
    cds_annotated_tx = parse_gtf(gtf_file, meta_d)    # CDS featureがすでにある転写産物は1, ないなら0

    now = datetime.now()
    timestamp = now.strftime("%y%m%d_%H%M%S")
    transcript_annotation_file = workdir + "/" + "transcript_annotation_" + timestamp + ".tsv"
    with open(transcript_annotation_file, mode="w") as f:
        for tx in transdecoder_result:
            # tx が入力GTFになければ飛ばす
            if not tx in cds_annotated_tx:
                continue
            # non-coding 予測の場合
            if tx.upper() in cpat_result:
                continue
            # すでにCDS feature が存在する転写産物の場合
            if cds_annotated_tx[tx] == 1:
                continue

            # 「coding potential 予測」かつ「CDS feature 無」の転写産物のみ処理
            start = transdecoder_result[tx]["cds_start"]
            end = transdecoder_result[tx]["cds_end"]
            coordinate = "nucleotide"
            feature = meta_d["cds"]
            source = "CPAT_TransDecoder"
            attributes = "."
            new_line = "\t".join([tx, str(start), str(end), coordinate, feature, source, attributes]) + "\n"
            f.write(new_line)

    # isoespy_makegtf.py を実行
    if gtf_file.endswith(".gtf"):
        cds_filled_gtf_file = workdir + "/" + os.path.basename(gtf_file)[:-4] + ".ORFpred.gtf"
    else:
        cds_filled_gtf_file = workdir + "/" + os.path.basename(gtf_file) + ".ORFpred.gtf"
    isoespy_makegtf(transcript_annotation_file, gtf_file, meta_file, cds_filled_gtf_file)


def main(args=None):
    parser = argparse.ArgumentParser(description="ORF prediction for unannotated transcripts by CPAT and TransDecoder.")
    parser.add_argument("--workdir", default=".", help="Working directory for isoespy_orfpred")
    parser.add_argument("--hexamer", required=True, help="Path to hexamer table (e.g., Human_Hexamer.tsv)")
    parser.add_argument("--model", required=True, help="Path to model data (e.g., Human_logitModel.RData)")
    parser.add_argument("--fasta", required=True, help="Input FASTA file of transcript sequences")
    parser.add_argument("--cpatoutprefix", required=True, help="Output prefix for CPAT results")
    parser.add_argument("--transdecoder_path", default="../TransDecoder-TransDecoder-v5.7.1", help="Path to TransDecoder executables (optional)")
    parser.add_argument("--gtf", required=True, help="GTF file INPUT")
    parser.add_argument("--meta", required=True, help="metadata")
    
    args = parser.parse_args()

    # Run CPAT
    run_cpat(args.workdir, args.hexamer, args.model, args.fasta, args.cpatoutprefix)
    print(">>> Finished CPAT.")

    # Run TransDecoder
    run_transdecoder(args.transdecoder_path, args.fasta, args.workdir)
    print(">>> Finished TransDecoder.")

    # ORF prediction
    transdecoder_outfile = args.workdir + "/" + os.path.basename(args.fasta) + ".transdecoder.bed"
    cpat_outfile = args.workdir + "/" + os.path.basename(args.cpatoutprefix) + ".no_ORF.txt"
    add_cds_to_gtf_from_cpat_transdecoder(args.workdir, transdecoder_outfile, cpat_outfile, args.gtf, args.meta)
    print(">>> Finished making complemented new GTF.")

if __name__ == "__main__":
    main()

