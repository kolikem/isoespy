library(DRIMSeq)
library(stageR)
library(yaml)

# ---- 引数の受け取り -------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
fileA <- args[1]  # count data (DRIMSeq input)
fileB <- args[2]  # samples data (DRIMSeq sample table)
fileC <- args[3]  # config (YAML: group_column, reference_group, coef_name)
fileD <- args[4]  # output p-values (stageR-adjusted)
fileE <- args[5]  # output proportions per sample

# ---- データ読み込み -------------------------------------------------
counts  <- read.delim(fileA, check.names = FALSE)
samples_df <- read.delim(fileB, check.names = FALSE)

config <- yaml::read_yaml(fileC)
group_col <- config$group_column         # 例: "group"
ref_group  <- config$reference_group     # 例: "Nontumor"
coef_name  <- config$coef_name           # 例: "groupHCC" みたいな係数名

# DRIMSeqオブジェクト作成
d <- dmDSdata(counts = counts, samples = samples_df)

# ---- サンプル情報の因子化とリファレンス設定 ------------------------
# samples(d) が使えないので d@samples を直接操作する
smp <- d@samples

# group列をfactor化してref_groupを基準にする
smp[[group_col]] <- factor(smp[[group_col]])
smp[[group_col]] <- stats::relevel(smp[[group_col]], ref = ref_group)

# 変更をオブジェクトに戻す
d@samples <- smp

# ---- フィルタリング -------------------------------------------------
# DRIMSeqのdmFilterには min_samps_gene_expr, min_samps_feature_expr などが必要
group_vec <- d@samples[[group_col]]
tab_group <- table(group_vec)

min_size <- min(tab_group)         # 各群の最小サンプル数
total_samples <- length(group_vec) # 全サンプル数

d <- dmFilter(
  d,
  min_samps_gene_expr    = total_samples*0.5,
  min_samps_feature_expr = min_size*0.5,
  min_gene_expr          = 10,
  min_feature_expr       = 10
)

# ---- Precision estimation ------------------------------------------
# デザイン行列を作る
# model.matrix(~ group, ...) を一般化して、group_col を使う
design_formula <- as.formula(paste("~", group_col))
design_full <- model.matrix(design_formula, data = d@samples)

set.seed(123)
d <- dmPrecision(d, design = design_full)

# ---- Proportion estimation (fit) -----------------------------------
d <- dmFit(d, design = design_full, verbose = 1)

# ---- 検定 (DIU test) ----------------------------------------------
# coef_name は config から与える（例: "group_colHCC" など）
d <- dmTest(d, coef = coef_name, verbose = 1)

# 結果をまとめる
res_gene <- DRIMSeq::results(d)                      # gene-level
res_gene <- res_gene[order(res_gene$pvalue), ]

res_tx   <- DRIMSeq::results(d, level = "feature")   # transcript-level

# ① NA をここで落としておく
res_gene_noNA <- res_gene[!is.na(res_gene$pvalue), ]
res_tx_noNA   <- res_tx[!is.na(res_tx$pvalue), ]

# もし全部NAならここで落とす
if (nrow(res_gene_noNA) == 0 || nrow(res_tx_noNA) == 0) {
  stop("All p-values are NA after filtering. Check filtering / design.")
}

# ---- stageR で multiple testing correction -------------------------
# ② stageR 用のオブジェクトをこの NA除去済みで作る
pScreen <- res_gene_noNA$pvalue
names(pScreen) <- res_gene_noNA$gene_id

pConfirmation <- matrix(res_tx_noNA$pvalue, ncol = 1)
rownames(pConfirmation) <- res_tx_noNA$feature_id

tx2gene <- res_tx_noNA[, c("feature_id", "gene_id")]

# ③ stageR 実行
stageRObj <- stageRTx(
  pScreen         = pScreen,
  pConfirmation   = pConfirmation,
  pScreenAdjusted = FALSE,
  tx2gene         = tx2gene
)

stageRObj <- stageWiseAdjustment(
  object  = stageRObj,
  method  = "dtu",
  alpha   = 0.05,
  allowNA = TRUE   # これはつけておいてOK
)

# ④ 調整済みp値を取得
padj <- getAdjustedPValues(
  stageRObj,
  order = TRUE,
  onlySignificantGenes = FALSE
)

#pScreen <- res_gene$pvalue
#names(pScreen) <- res_gene$gene_id
#
#pConfirmation <- matrix(res_tx$pvalue, ncol = 1)
#rownames(pConfirmation) <- res_tx$feature_id
#
#tx2gene <- res_tx[, c("feature_id", "gene_id")]
#
#stageRObj <- stageRTx(
#  pScreen = pScreen,
#  pConfirmation = pConfirmation,
#  pScreenAdjusted = FALSE,
#  tx2gene = tx2gene
#)
#
#stageRObj <- stageWiseAdjustment(
#  object = stageRObj,
#  method = "dtu",
#  alpha = 0.05
#)
#
#padj <- getAdjustedPValues(
#  stageRObj,
#  order = TRUE,
#  onlySignificantGenes = FALSE
#)

# ---- 出力1: p値/調整済みp値 ----------------------------------------
write.table(
  padj,
  file = fileD,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---- 出力2: 各サンプルにおけるアイソフォーム使用率 -------------------
# proportions(d) は isoform usage (各sampleでtranscriptが占める割合)
# これはDRIMSeq側の関数でOK
usage_mat <- DRIMSeq::proportions(d)

write.table(
  usage_mat,
  file = fileE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

