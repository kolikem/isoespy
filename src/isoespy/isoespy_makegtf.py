import subprocess
import argparse
import os
try:
    import proteinToGenome
except ModuleNotFoundError:
    import isoespy.proteinToGenome

### GTF
#	chr6	transdecoder	CDS	47712501	47712608	0.000000	+	.	gene_id "-"; transcript_id "ADGRF4_novel_tx_1"; gene_symbol "ADGRF4"; SPLICE_CATEGORY "novel_exon_combination";


### GTF prep format (transcript-based annotation format)
#	#Transcript		from	to	coordinate	feature			source		attributes
#	ADGRF4_novel_tx_1	252	500	amino_acid	protein_domain		pfam		domain_id "PF00002.29"; domain_name "7tm_2"; evalue "3.2e-33";
#	ADGRF4_novel_tx_1	199	242	amino_acid	protein_domain		pfam		domain_id "PF01825.27"; domain_name "GPS"; evalue; "4e-10";
#	ENST00000283303.3	252	500	amino_acid	protein_domain		pfam		domain_id "PF00002.29"; domain_name "7tm_2"; evalue "5.5e-33";
#	ENST00000283303.3	347	390	amino_acid	protein_domain		pfam		domain_id "PF01825.27"; domain_name "GPS"; evalue "5.4e-10";


def parse_meta(meta_file):
    ### 1. GTF meta
    # 初期値
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

    return transcript_id, exon, cds


def parse_annotfile(annot_file):
    annot_data = {}
    feat_nbr = 0
    with open(annot_file) as f:
        for line in f:
            if line.startswith("#"):   # ヘッダー行は"#"で開始する
                continue
            feat_nbr += 1
            cols = line.strip().split('\t')

            tx_name = cols[0]
            featfrom = int(cols[1])
            featto = int(cols[2])
            coord = cols[3]
            feature = cols[4]
            source = cols[5]
            attribute = cols[6]
            if attribute == ".":
                attribute = ""

            # 処理
            if not tx_name in annot_data:
                annot_data[tx_name] = []
            annot_data[tx_name].append({"coordinate": coord, "from": featfrom, "to": featto, "feature": feature, "source": source, "attribute": attribute, "feature_id": feat_nbr})
    return annot_data


def parse_gtf(gtf_file, metadata_file, TRANSCRIPT_ID, EXON, CDS):
    """GTFファイルをパースし、転写産物ごとのcds情報を取得"""
    # TRANSCRIPT_ID : GTFファイル内におけるattributes 列の転写産物名フィールドのキー名
    # CDS           : GTFファイル内におけるfeature 列のCDS featureの名前
    # EXON          : GTFファイル内におけるfeature 列のexon featureの名前

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
            transcript_id = attr_dict.get(TRANSCRIPT_ID)

            # 処理
            if transcript_id:
                if feature == CDS:
                    if transcript_id not in gtf_data:
                        gtf_data[transcript_id] = {"chrom": chrom, "strand": strand, "exon": [], "cds": [], "score": score, "strand": strand, "frame": frame, "attributes": attributes}
                    gtf_data[transcript_id]["cds"].append((start, end))
                if feature == EXON:
                    if transcript_id not in gtf_data:
                        gtf_data[transcript_id] = {"chrom": chrom, "strand": strand, "exon": [], "cds": [], "score": score, "strand": strand, "frame": frame, "attributes": attributes}
                    gtf_data[transcript_id]["exon"].append((start, end))

    # exon を座標順に並び替え
    for tx in gtf_data:
        gtf_data[tx]["exon"] = sorted(gtf_data[tx]["exon"], key=lambda x: x[0])
    # CDS を座標順に並び替え
    for tx in gtf_data:
        gtf_data[tx]["cds"] = sorted(gtf_data[tx]["cds"], key=lambda x: x[0])

    # 塩基配列ベースのexon総延長(転写産物長)を計算する
    totallen_exon = {}
    for tx in gtf_data:
        exon_len = 0
        for i in gtf_data[tx]["exon"]:
            exon_len += i[1] - i[0] + 1
            totallen_exon[tx] = exon_len
    # 塩基配列ベースのCDS総延長を計算する
    totallen_cds = {}
    for tx in gtf_data:
        cds_len = 0
        for i in gtf_data[tx]["cds"]:
            cds_len += i[1] - i[0] + 1
            totallen_cds[tx] = cds_len
    for tx in gtf_data:
        if not tx in totallen_cds:
            totallen_cds[tx] = 0

    return gtf_data, totallen_exon, totallen_cds


def parse_frame(frame_file):
    frame_d = {}
    if frame_file == None:
        return {}

    else:
        with open(frame_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                frame_d[cols[0]] = int(cols[1])
        return frame_d


def ConvertAnnotationDataToGTF(annot_data, total_exon_d, total_cds_d, gtf_data, outfile_name, frame_d):
    for tx in annot_data:
        if tx in gtf_data:
            coord_gen_cds = gtf_data[tx]["cds"]
            coord_gen_exon = gtf_data[tx]["exon"]
            strand_gen = gtf_data[tx]["strand"]
            chrom = gtf_data[tx]["chrom"]
            score = gtf_data[tx]["score"]
            total_exon = total_exon_d[tx]
            total_cds = total_cds_d[tx]

            # for each feature for the transcript
            for feature_i in annot_data[tx]:
                coord_i = feature_i["coordinate"]
                from_i =  feature_i["from"]
                to_i = feature_i["to"]
                source = feature_i["source"]
                feature = feature_i["feature"]
                attributes_annot = feature_i["attribute"]
                feature_id = str(feature_i["feature_id"])

                # backmap
                if coord_i == "nucleotide":
                    backcoord_i = proteinToGenome.backmap_nuc(coord_gen_exon, tx, from_i, to_i, total_exon, strand_gen)		# ゲノム上にマップした塩基配列範囲
                elif coord_i == "amino_acid":
                    backcoord_i = proteinToGenome.backmap_aa(coord_gen_cds, tx, from_i, to_i, total_cds, strand_gen, frame_d)	# ゲノム上にマップしたアミノ酸配列範囲

                # GTFにfeature_iを追加
                for j in backcoord_i:
                    start, end = str(j[0]), str(j[1])
                    frame = "."

                    attributes_gtf = gtf_data[tx]["attributes"]                    
                    # attributeフィールドを統合 (annotation --> GTF)
                    attr_dict_annot = {k.strip(): v.strip('"') for k, v in (item.split() for item in attributes_annot.split(';') if item)}
                    attr_dict_gtf = {k.strip(): v.strip('"') for k, v in (item.split() for item in attributes_gtf.split(';') if item)}
                    for a in attr_dict_annot:
                        if not a in attr_dict_gtf:
                            # GTFのattributeフィールドに元々存在しないannotationファイル内のattributeを加える
                            attr_dict_gtf[a] = attr_dict_annot[a]
                    attributes_l = []
                    for a in attr_dict_gtf:
                        attributes_l.append(a+' '+'"'+attr_dict_gtf[a]+'"'+';')
                    attributes = ' '.join(attributes_l)
                    # 最後にfeature_id を追記
                    attributes += ' feature_id '+'"'+feature_id+'"'+';'
                    
                    with open(outfile_name, mode="a") as f:
                        f.write("\t".join([chrom, source, feature, start, end, score, strand_gen, frame, attributes]))
                        f.write("\n")


def outfileName(annotation_file, gtf_file):
    annot_base = os.path.basename(annotation_file)
    gtf_base = os.path.basename(gtf_file)
    annot_dir = os.path.dirname(annotation_file)
    return annot_dir + "/" + gtf_base + "_" + annot_base + "_" + "isoespy_makegtf.gtf"


def outfileCopy(gtf_file, outfile_name):
    subprocess.run(["rm", "-f", outfile_name])
    subprocess.run(["cp", gtf_file, outfile_name])


def separate_annotData(annot_data):
    """ GTF prep データをCDSの行とそれ以外のアノテーションの部分に分ける.
        CDSは先にGTF化する.
    """
    annot_data_cds = {}
    annot_data_else = {}
    for tx, features in annot_data.items():
        cds_list = [f for f in features if f.get("feature") == "CDS"]
        else_list = [f for f in features if f.get("feature") != "CDS"]
        if cds_list:
            annot_data_cds[tx] = cds_list
        if else_list:
            annot_data_else[tx] = else_list
    return annot_data_cds, annot_data_else


def isoespy_makegtf(gtfprep_file, gtf_file, meta_file, outfile_name, frame_file):
    # Metadata
    TRANSCRIPT_ID, EXON, CDS = parse_meta(meta_file)

    # frame_data
    frame_d = parse_frame(frame_file)

    # gtfprep data の分割
    annot_data = parse_annotfile(gtfprep_file)
    annot_data_cds, annot_data_else = separate_annotData(annot_data)

    # CDS --> GTF
    gtf_data, total_exon_d, total_cds_d = parse_gtf(gtf_file, meta_file, TRANSCRIPT_ID, EXON, CDS)
    # Get output file name
    if outfile_name == None:
        outfile_name = outfileName(gtfprep_file, gtf_file)
    # Initiate output file
    outfileCopy(gtf_file, outfile_name)
    # Make GTF
    ConvertAnnotationDataToGTF(annot_data_cds, total_exon_d, total_cds_d, gtf_data, outfile_name, frame_d)

    # non-CDS feature(s) --> GTF
    gtf_data, total_exon_d, total_cds_d = parse_gtf(outfile_name, meta_file, TRANSCRIPT_ID, EXON, CDS)
    # Make GTF
    ConvertAnnotationDataToGTF(annot_data_else, total_exon_d, total_cds_d, gtf_data, outfile_name, frame_d)


def main(args=None):
    # parse
    parser = argparse.ArgumentParser(description='isoespy_makegtf()')
    parser.add_argument('-p', '--gtfprep_file', required=True, type=str, default=None, help='GTF prep formatted file as input')
    parser.add_argument('-g', '--gtf_file', required=True, type=str, default=None, help='GTF file as input')
    parser.add_argument('-m', '--meta_file', required=True, type=str, default=None, help='meta file as input')
    parser.add_argument('-o', '--output', required=True, type=str, default=None, help='output file')
    parser.add_argument('-f', '--frame_file', type=str, default=None, help='Transcript translation frame file. Required when coordinate in GTF prep file contains amino_acid whose frame of coordinate may not be 0.')

    args = parser.parse_args()

    isoespy_makegtf(args.gtfprep_file, args.gtf_file, args.meta_file, args.output, args.frame_file)


if __name__ == "__main__":
    main()

