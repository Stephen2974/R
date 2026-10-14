library(Seurat)
library(data.table)
library(R.utils)  # fread 直接读 .gz 需要
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(scDblFinder)
library(SingleCellExperiment)
library(Matrix)
library(harmony)

#####1. 读取目录下 CSV.gz 表达矩阵（基因 x 细胞）####
# 数据在 G:/R，无需先解压 .gz
setwd("G:/R")

# 获取所有样本文件（如 GSM5573466_sample1.csv.gz）
files <- list.files("G:/R", pattern = "\\.csv\\.gz$", full.names = TRUE)
files <- sort(files)

# 创建一个空的列表
seurat_list <- list()

# 读取数据并创建Seurat对象
for (f in files) {
  
  # 样本名保留完整文件名（去后缀）：GSM5573466_sample1.csv.gz -> GSM5573466_sample1
  sample <- sub("\\.csv\\.gz$", "", basename(f))
  # 单独抽出 GSM 登录号，便于对照 GEO
  gsm_id <- sub("_.*$", "", sample)
  
  message("正在读取: ", basename(f), " -> ", sample)
  
  # 读入压缩 CSV：第 1 列基因名，其余列为细胞
  dt <- fread(f)
  genes <- as.character(dt[[1]])
  counts <- as.matrix(dt[, -1, with = FALSE])
  rownames(counts) <- make.unique(genes)
  counts <- Matrix(counts, sparse = TRUE)
  rm(dt)
  
  # 创建Seurat对象（orig.ident = 完整样本名）
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample,
    min.features = 200,
    min.cells = 3
  )
  
  # 元数据：文件名能确定的字段；Patient/Group 等有表后再补
  seurat_obj$Sample <- sample
  seurat_obj$GSM <- gsm_id
  
  # 添加到列表中
  seurat_list[[sample]] <- seurat_obj
  rm(counts, seurat_obj)
}

# 合并所有对象
if (length(seurat_list) == 1) {
  seu <- seurat_list[[1]]
} else {
  seu <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list)
  )
}
# 合并刚结束、质控还没开始，此时每个细胞都在
# 先补上病人、部位和分组，后面的质控和子集会一直带着这三列
sample_meta <- data.frame(
  # sample1 到 sample40。下面三列必须按这个顺序写，第 1 个值对应 sample1
  sample_no = paste0("sample", 1:40),
  Patient = c(  # 病人编号。同一个人的肿瘤和正常用同一个编号
    "NGCII518","NGCII518","NGCII519","NGCII520","NGCII520",
    "NGCII521","NGCII521","NGCII524","NGCII525","NGCII525",
    "NGCII527","NGCII527","NGCII529","NGCII502","NGCII511",
    "NGCII499","NGCII509","NGCII498","NGCII011","NGCII015",
    "NGCII513","NGCII513","NGCII514","NGCII514","NGCII512",
    "NGCII512","NGCII510","NGCII531","NGCII522","NGCII533",
    "NGCII536","NGCII539","NGCII540","NGCII541","NGCII538",
    "NGCII538","NGCII538","NGCII538","NGCII543","NGCII545"
  ),
  Tissue = c(   # 取样部位：原发胃，或腹膜
    "原发胃","原发胃","原发胃","原发胃","原发胃",
    "原发胃","原发胃","原发胃","原发胃","原发胃",
    "原发胃","原发胃","原发胃","原发胃","原发胃",
    "原发胃","原发胃","原发胃","腹膜","腹膜",
    "原发胃","原发胃","原发胃","原发胃","原发胃",
    "原发胃","原发胃","原发胃","原发胃","原发胃",
    "原发胃","原发胃","原发胃","原发胃","原发胃",
    "原发胃","腹膜","腹膜","原发胃","原发胃"
  ),
  Group = c(    # 肿瘤，或配对的正常胃组织
    "正常","肿瘤","肿瘤","正常","肿瘤",
    "正常","肿瘤","肿瘤","正常","肿瘤",
    "正常","肿瘤","肿瘤","肿瘤","肿瘤",
    "肿瘤","肿瘤","肿瘤","肿瘤","肿瘤",
    "正常","肿瘤","正常","肿瘤","正常",
    "肿瘤","肿瘤","肿瘤","肿瘤","肿瘤",
    "正常","肿瘤","肿瘤","肿瘤","正常",
    "肿瘤","正常","肿瘤","肿瘤","肿瘤"
  ),
  stringsAsFactors = FALSE
)

# Sample 现在是 GSM5573466_sample1。去掉前面的 GSM 编号，只留 sample1，才能和对照表对齐
seu$sample_no <- sub("^.*_", "", seu$Sample)
# 每个细胞对应对照表的哪一行。对不上时是 NA
hit <- match(seu$sample_no, sample_meta$sample_no)
stopifnot(!any(is.na(hit)))   # 有 NA 就停下来，说明有样本没写进对照表
# 按行号把病人、部位、分组写到每个细胞上
seu$Patient <- sample_meta$Patient[hit]
seu$Tissue  <- sample_meta$Tissue[hit]
seu$Group   <- sample_meta$Group[hit]

# 核对：每个样本只能是肿瘤或正常；腹膜应只有 sample19、20、38 的肿瘤和 sample37 的正常
table(seu$sample_no, seu$Group)
table(seu$Tissue, seu$Group)
# 验证结果
head(seu@meta.data)
table(seu@meta.data$orig.ident)

#####2. 质控（全部样本在 seu 上一次完成）####
# 线粒体基因比例
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# 红细胞比例
HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
HB.genes.use <- CaseMatch(HB.genes, rownames(seu))
seu[["percent.HB"]] <- PercentageFeatureSet(seu, features = HB.genes.use)

# 查看相关性
FeatureScatter(seu, "nCount_RNA", "percent.mt", group.by = "orig.ident")
FeatureScatter(seu, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")

# 质控指标小提琴图（质控前）
theme.set2 <- theme(axis.title.x = element_blank())
plot.features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group <- "orig.ident"
plots <- list()
for (i in seq_along(plot.features)) {
  plots[[i]] <- VlnPlot(seu, group.by = group, pt.size = 0,
                        features = plot.features[i]) + theme.set2 + NoLegend()
}
violin_before <- wrap_plots(plots = plots, nrow = 2)
violin_before
ggsave("violin_seu_before_qc.pdf", plot = violin_before, width = 20, height = 8)
dim(seu)

# 分位数，便于根据本数据集调整阈值
quantile(seu$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(seu$nFeature_RNA, seq(0.9, 1, 0.01))
quantile(seu$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(seu$nCount_RNA, seq(0.9, 1, 0.01))
quantile(seu$percent.mt, seq(0.9, 1, 0.01))
quantile(seu$percent.HB, seq(0.9, 1, 0.01))

# 质控阈值（可先看小提琴图/分位数后再改；参考原 GEO 处理约 500–6000、mt<20）
minGene <- 500
maxGene <- 6000
minUMI <- 600
pctMT <- 20
pctHB <- 1

seu <- subset(
  seu,
  subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
    nCount_RNA > minUMI & percent.mt < pctMT & percent.HB < pctHB
)

# 质控后小提琴图
plots <- list()
for (i in seq_along(plot.features)) {
  plots[[i]] <- VlnPlot(seu, group.by = group, pt.size = 0,
                        features = plot.features[i]) + theme.set2 + NoLegend()
}
violin_after <- wrap_plots(plots = plots, nrow = 2)
violin_after
ggsave("violin_seu_after_qc.pdf", plot = violin_after, width = 20, height = 8)
dim(seu)

#####3. 去双细胞（按样本估计）####
# Seurat v5：合并后先 JoinLayers
seu <- JoinLayers(seu)

# 转成 SingleCellExperiment，按样本跑 scDblFinder
sce <- as.SingleCellExperiment(seu)
sce <- scDblFinder(sce, samples = sce$orig.ident)
table(sce$scDblFinder.class)

# 转回 Seurat 并只保留 singlet
seu <- as.Seurat(sce, counts = "counts", data = NULL)
table(seu$scDblFinder.class)
seu_clean <- subset(seu, subset = scDblFinder.class == "singlet")
message("去双细胞后保留细胞数: ", ncol(seu_clean))

seu_clean$Is_Double <- seu_clean$scDblFinder.class
seu_clean$Double_Score <- seu_clean$scDblFinder.score
head(seu_clean@meta.data[, c("Double_Score", "Is_Double")])

#####4. 归一化、细胞周期、降维整合####
seu_clean <- NormalizeData(seu_clean)
seu_clean <- FindVariableFeatures(seu_clean, selection.method = "vst", nfeatures = 2000)

seu_clean <- CellCycleScoring(
  seu_clean,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = FALSE
)

# 先不回归细胞周期。S.Score 和 G2M.Score 先留在元数据里，后面看它们是否真的主导分群
seu_clean <- ScaleData(seu_clean)
saveRDS(seu_clean, file = "G:/R/seu_clean_quxibaozhouqi.rds")

seu_clean <- RunPCA(seu_clean)

# 此时 PCA 还没去掉周期。看前 10 个 PC 和 S / G2M 的相关
# 相关弱就用现在这份结果；若某个 PC 和周期分数明显一起变，再单独做一版回归
pca_emb <- Embeddings(seu_clean, reduction = "pca")
cor_s <- cor(pca_emb[, 1:10], seu_clean$S.Score, method = "spearman")
cor_g2m <- cor(pca_emb[, 1:10], seu_clean$G2M.Score, method = "spearman")
cor_df <- data.frame(
  PC = paste0("PC", 1:10),
  cor_S = round(as.numeric(cor_s), 3),
  cor_G2M = round(as.numeric(cor_g2m), 3)
)
cor_df

# Harmony 按样本去批次（harmony 2.x 用 reduction.use，不要写 reduction）
seu_clean <- RunHarmony(
  seu_clean,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  dims.use = 1:30
)

# 聚类 + UMAP
seu_clean <- FindNeighbors(seu_clean, dims = 1:30, reduction = "harmony")
seu_clean <- FindClusters(seu_clean, resolution = 0.3)
seu_clean <- RunUMAP(seu_clean, dims = 1:30, reduction = "harmony")

# 按肿瘤/正常、原发胃/腹膜、病人着色。前两张保留图例；病人太多，第三张不画图例
ggsave("umap_group.pdf",
       DimPlot(seu_clean, reduction = "umap", group.by = "Group"),
       width = 8, height = 6)
ggsave("umap_tissue.pdf",
       DimPlot(seu_clean, reduction = "umap", group.by = "Tissue"),
       width = 8, height = 6)
ggsave("umap_patient.pdf",
       DimPlot(seu_clean, reduction = "umap", group.by = "Patient") + NoLegend(),
       width = 8, height = 6)
# 簇编号图，后面做大类注释时对照这张
ggsave("umap_cluster.pdf",
       DimPlot(seu_clean, reduction = "umap", label = TRUE),
       width = 8, height = 6)
saveRDS(seu_clean, file = "G:/R/seu_clean_harmony_umap.rds")
