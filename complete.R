library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(scDblFinder)
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(Matrix)
library(harmony)

#####1. 读取目录下 CSV.gz 表达矩阵（基因 x 细胞）####
# 数据在 D:/R，无需先解压 .gz
setwd("D:/R")

# 获取所有样本文件（如 GSM5573466_sample1.csv.gz）
files <- list.files("D:/R", pattern = "\\.csv\\.gz$", full.names = TRUE)
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

seu_clean <- ScaleData(seu_clean, vars.to.regress = c("S.Score", "G2M.Score"))
saveRDS(seu_clean, file = "D:/R/seu_clean_quxibaozhouqi.rds")

seu_clean <- RunPCA(seu_clean)

# 检查前 10 个 PC 与细胞周期分数的相关
pca_emb <- Embeddings(seu_clean, reduction = "pca")
cor_s <- cor(pca_emb[, 1:10], seu_clean$S.Score, method = "spearman")
cor_g2m <- cor(pca_emb[, 1:10], seu_clean$G2M.Score, method = "spearman")
cor_df <- data.frame(
  PC = paste0("PC", 1:10),
  cor_S = round(as.numeric(cor_s), 3),
  cor_G2M = round(as.numeric(cor_g2m), 3)
)
cor_df

# Harmony 按样本去批次
seu_clean <- RunHarmony(
  seu_clean,
  group.by.vars = "orig.ident",
  reduction = "pca",
  dims.use = 1:30
)

# 聚类 + UMAP
seu_clean <- FindNeighbors(seu_clean, dims = 1:30, reduction = "harmony")
seu_clean <- FindClusters(seu_clean, resolution = 0.3)
seu_clean <- RunUMAP(seu_clean, dims = 1:30, reduction = "harmony")

# 简单查看
DimPlot(seu_clean, reduction = "umap", group.by = "orig.ident", label = FALSE) + NoLegend()
DimPlot(seu_clean, reduction = "umap", label = TRUE)
saveRDS(seu_clean, file = "D:/R/seu_clean_harmony_umap.rds")
