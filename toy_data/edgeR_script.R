library(edgeR)

# 引数
args <- commandArgs(trailingOnly=TRUE)
fileA <- readLines(args[1])
fileB <- args[2]
fileC <- args[3]
fileD <- args[4]

# [sample] セクションの開始位置を探す
sample_start <- which(fileA == "[sample]") + 1

# 次のブロックの開始位置を探す（例：[config] や他のブロック）
next_block_start <- which(grepl("^\\[.*\\]$", fileA))

# [sample] セクションの終了位置を決定
if (length(next_block_start) > 0) {
  # [sample] セクション以降にブロックがある場合
  next_block_start <- next_block_start[next_block_start > sample_start]  # [sample] より後のブロックを探す
  if (length(next_block_start) > 0) {
    sample_end <- min(next_block_start) - 1
  } else {
    sample_end <- length(fileA)
  }
} else {
  # 他のブロックがない場合は、ファイルの最後まで
  sample_end <- length(fileA)
}

# [sample] セクションのみを抽出
sample_lines <- fileA[sample_start:sample_end]

# controlとtargetのパース
levels <- fileA[sample_start]
split_levels <- strsplit(levels, ",")[[1]]
parsed_values <- setNames(sapply(strsplit(split_levels, "="), `[`, 2),
                          sapply(strsplit(split_levels, "="), `[`, 1))
control_group <- parsed_values["#control"]	# N
target_group <- parsed_values["target"]		# T

# データフレームに変換
sample_lines <- sample_lines[!grepl("^#", sample_lines)]
metadata <- read.table(text = paste(sample_lines, collapse = "\n"), header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("sample", "group"))

# ファイルB（カウントデータ）の読み込み
counts <- read.table(fileB, header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1, check.names=FALSE)

# サンプル名を取得
sample_names <- colnames(counts)
#print(sample_names)

# メタデータの群情報を取得
group <- metadata$group[match(sample_names, metadata$sample)]	# sample_namesとgroupはカウントデータの順序
group <- factor(group, levels=c(control_group, target_group))

# 遺伝子発現データをDGElistオブジェクトに変換
DGE <- DGEList(counts=counts, group=group)

# TMM正規化係数
DGE <- calcNormFactors(DGE)

# edgeR CPM
cpm_values <- cpm(DGE, normalized.lib.sizes=TRUE)

# statistical test
design <- model.matrix(~group)
DGE <- estimateDisp(DGE, design)
fit_DGE <- glmQLFit(DGE, design)
qlf_DGE <- glmQLFTest(fit_DGE)

# result (tests)
result <- topTags(qlf_DGE, n=Inf)
result_df <- result$table
result_df$gene <- rownames(result_df)
write.csv(result_df, fileC, quote=FALSE, row.names=FALSE)

# result (CPM)
write.table(t(c("gene", as.character(colnames(cpm_values)))), file=fileD, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(cpm_values, file=fileD, sep="\t", quote=FALSE, row.names=TRUE, col.names=FALSE, append=TRUE)

