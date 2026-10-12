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
#####1. 读取目录下 CSV.gz 表达矩阵（基因 x 细胞）####
# 数据在 D:/R，无需先解压 .gz
library(Matrix)
setwd("D:/R")

# 获取所有样本文件（如 GSM5573466_sample1.csv.gz）
files <- list.files("D:/R", pattern = "\\.csv\\.gz$", full.names = TRUE)
files <- sort(files)

# 创建一个空的列表
seurat_list <- list()

# 读取数据并创建Seurat对象
for (f in files) {
  
  # 样本名：GSM5573466_sample1.csv.gz -> sample1
  sample <- sub("\\.csv\\.gz$", "", basename(f))
  sample <- sub("^GSM\\d+_", "", sample)
  
  message("正在读取: ", basename(f), " -> ", sample)
  
  # 读入压缩 CSV：第 1 列基因名，其余列为细胞
  dt <- fread(f)
  genes <- as.character(dt[[1]])
  counts <- as.matrix(dt[, -1, with = FALSE])
  rownames(counts) <- make.unique(genes)
  counts <- Matrix(counts, sparse = TRUE)
  rm(dt)
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample,
    min.features = 200,
    min.cells = 3
  )
  
  # 写入样本信息（原 Patient/Day 命名不适用当前 GEO 文件名，后续有 metadata 再补）
  seurat_obj$Sample <- sample
  
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

# ------------------------------------------------------------
# 说明：读入已改为单列表 + 合并对象 seu。
# 下面质控仍按旧的 HC / PT283 / PT406 三组写，下一步再改分组逻辑。
# ------------------------------------------------------------

#质量控制
#线粒体基因比例
HC[["percent.mt"]] <- PercentageFeatureSet(HC, pattern = "^MT-")
PT283[["percent.mt"]] <- PercentageFeatureSet(PT283, pattern = "^MT-")
PT406[["percent.mt"]] <- PercentageFeatureSet(PT406, pattern = "^MT-")
#红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes1 <- CaseMatch(HB.genes, rownames(HC))
HB.genes2 <- CaseMatch(HB.genes, rownames(PT283))
HB.genes3 <- CaseMatch(HB.genes, rownames(PT406))
HC[["percent.HB"]]<-PercentageFeatureSet(HC, features=HB.genes1) 
PT283[["percent.HB"]]<-PercentageFeatureSet(PT283, features=HB.genes2) 
PT406[["percent.HB"]]<-PercentageFeatureSet(PT406, features=HB.genes3) 
#查看相关性
FeatureScatter(HC, "nCount_RNA", "percent.mt", group.by = "orig.ident")
FeatureScatter(PT283, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")
FeatureScatter(PT406, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")
#查看质控指标
#设置绘图元素
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group = "orig.ident"
#质控前小提琴图
plots = list()
for(i in c(1:length(plot.featrures))){
  plots[[i]] = VlnPlot(HC, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin_HC <- wrap_plots(plots = plots, nrow=2)  
violin_HC
for(i in c(1:length(plot.featrures))){
  plots[[i]] = VlnPlot(PT283, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin_PT283 <- wrap_plots(plots = plots, nrow=2)  
violin_PT283
for(i in c(1:length(plot.featrures))){
  plots[[i]] = VlnPlot(PT406, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin_PT406 <- wrap_plots(plots = plots, nrow=2)  
violin_PT406
#保存
ggsave("violin_HC_before1_qc.pdf", plot = violin_HC, width = 14, height = 8) 
ggsave("violin_PT283_before1_qc.pdf", plot = violin_PT283, width = 14, height = 8) 
ggsave("violin_PT406_before1_qc.pdf", plot = violin_PT406, width = 14, height = 8) 
dim(HC)
dim(PT283)
dim(PT406)
#设置质控指标
quantile(HC$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(HC$nFeature_RNA, seq(0.9, 1, 0.01))
quantile(PT283$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(PT283$nFeature_RNA, seq(0.9, 1, 0.01))
quantile(PT406$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(PT406$nFeature_RNA, seq(0.9, 1, 0.01))
#plots[[1]] + geom_hline(yintercept = 500) + geom_hline(yintercept = 4500)

quantile(HC$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(HC$nCount_RNA, seq(0.9, 1, 0.01))
quantile(PT283$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(PT283$nCount_RNA, seq(0.9, 1, 0.01))
quantile(PT406$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(PT406$nCount_RNA, seq(0.9, 1, 0.01))
#plots[[2]] + geom_hline(yintercept = 22000)
quantile(HC$percent.mt, seq(0.9, 1, 0.01))
quantile(PT283$percent.mt, seq(0.9, 1, 0.01))
quantile(PT406$percent.mt, seq(0.9, 1, 0.01))
#plots[[3]] + geom_hline(yintercept = 20)
quantile(HC$percent.HB, seq(0.9, 1, 0.01))
quantile(PT283$percent.HB, seq(0.9, 1, 0.01))
quantile(PT406$percent.HB, seq(0.9, 1, 0.01))
#plots[[4]] + geom_hline(yintercept = 1)

#设置质控标准
#基因
minGene_HC=500
maxGene_HC=3000
minGene_PT283=500
maxGene_PT283=3500
minGene_PT406=500
maxGene_PT406=3000
#counts
minUMI=600
#线粒体
pctMT=20
#血细胞
pctHB=1

#数据质控并绘制小提琴图
HC <- subset(HC, subset = nFeature_RNA > minGene_HC & nFeature_RNA < maxGene_HC &
               nCount_RNA > minUMI  & percent.mt < pctMT & percent.HB < pctHB)
PT283 <- subset(PT283, subset = nFeature_RNA > minGene_PT283 & nFeature_RNA < maxGene_PT283 &
                  nCount_RNA > minUMI &  percent.mt < pctMT & percent.HB < pctHB)
PT406 <- subset(PT406, subset = nFeature_RNA > minGene_PT406 & nFeature_RNA < maxGene_PT406 &
                  nCount_RNA > minUMI &  percent.mt < pctMT & percent.HB < pctHB)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(HC, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin_HC <- wrap_plots(plots = plots, nrow=2)    
violin_HC
dim(HC)
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(PT283, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin_PT283 <- wrap_plots(plots = plots, nrow=2)    
violin_PT283
dim(PT283)
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(PT406, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin_PT406 <- wrap_plots(plots = plots, nrow=2)    
violin_PT406
dim(PT406)
#HC
# 1. 合并 layers（仅 Seurat v5 需要）
HC <- JoinLayers(HC)
# 2. 转成 SingleCellExperiment
sce <- as.SingleCellExperiment(HC)
# 3. 运行 scDblFinder（单样本不需要 clusters 和 samples）
sce <- scDblFinder(sce)
# 4. 查看双细胞分类
table(sce$scDblFinder.class)
# 5. 转回 Seurat
HC <- as.Seurat(sce, counts = "counts", data = NULL)
# 6. 确认 meta.data 中有 scDblFinder.class
table(HC$scDblFinder.class)
# 7. 过滤双细胞
HC_clean <- subset(HC, subset = scDblFinder.class == "singlet")
print(paste("过滤后保留的细胞数:", ncol(HC_clean)))
#将列名改为"Double_Score"和"Is_Double"
colnames(HC_clean@meta.data)
HC_clean$Is_Double <- HC_clean$scDblFinder.class
HC_clean$Double_Score <- HC_clean$scDblFinder.score
#查看DoubletFinder分析结果
head(HC_clean@meta.data[, c("Double_Score", "Is_Double")])
#PT283
PT283 <- JoinLayers(PT283)
# 2. 转成 SingleCellExperiment
sce2 <- as.SingleCellExperiment(PT283)
# 3. 运行 scDblFinder（单样本不需要 clusters 和 samples）
sce2 <- scDblFinder(sce2)
# 4. 查看双细胞分类
table(sce2$scDblFinder.class)
# 5. 转回 Seurat
PT283 <- as.Seurat(sce2, counts = "counts", data = NULL)
# 6. 确认 meta.data 中有 scDblFinder.class
table(PT283$scDblFinder.class)
# 7. 过滤双细胞
PT283_clean <- subset(PT283, subset = scDblFinder.class == "singlet")
print(paste("过滤后保留的细胞数:", ncol(PT283_clean)))
#将列名改为"Double_Score"和"Is_Double"
colnames(PT283_clean@meta.data)
PT283_clean$Is_Double <- PT283_clean$scDblFinder.class
PT283_clean$Double_Score <- PT283_clean$scDblFinder.score
#查看DoubletFinder分析结果
head(PT283_clean@meta.data[, c("Double_Score", "Is_Double")])
#PT406
PT406 <- JoinLayers(PT406)
# 2. 转成 SingleCellExperiment
sce3 <- as.SingleCellExperiment(PT406)
# 3. 运行 scDblFinder（单样本不需要 clusters 和 samples）
sce3 <- scDblFinder(sce3)
# 4. 查看双细胞分类
table(sce3$scDblFinder.class)
# 5. 转回 Seurat
PT406 <- as.Seurat(sce3, counts = "counts", data = NULL)
# 6. 确认 meta.data 中有 scDblFinder.class
table(PT406$scDblFinder.class)
# 7. 过滤双细胞
PT406_clean <- subset(PT406, subset = scDblFinder.class == "singlet")
print(paste("过滤后保留的细胞数:", ncol(PT406_clean)))
#将列名改为"Double_Score"和"Is_Double"
colnames(PT406_clean@meta.data)
PT406_clean$Is_Double <- PT406_clean$scDblFinder.class
PT406_clean$Double_Score <- PT406_clean$scDblFinder.score
#查看DoubletFinder分析结果
head(PT406_clean@meta.data[, c("Double_Score", "Is_Double")])
#合并样本
combined <- merge(
  x = HC_clean,
  y = list(PT283_clean, PT406_clean),
  add.cell.ids = c("HC", "PT283", "PT406"),
  project = "liver_transplant"
)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
# 1. 细胞周期评分（必须放在 NormalizeData 之后）
combined <- CellCycleScoring(
  combined,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = FALSE
)

# 2. ScaleData 时回归掉 S.Score 和 G2M.Score
combined <- ScaleData(combined, vars.to.regress = c("S.Score", "G2M.Score"))
saveRDS(combined, file = "D:/liver2/combined_clean(quxibaozhouqi).rds")
# 3. 然后才是 PCA
combined <- RunPCA(combined)
# 提取 PCA 坐标
pca_emb <- Embeddings(combined, reduction = "pca")

# 计算前 10 个主成分与 S.Score、G2M.Score 的相关系数
cor_s   <- cor(pca_emb[, 1:10], combined$S.Score,   method = "spearman")
cor_g2m <- cor(pca_emb[, 1:10], combined$G2M.Score, method = "spearman")

# 汇总
cor_df <- data.frame(
  PC = paste0("PC", 1:10),
  cor_S   = round(as.numeric(cor_s), 3),
  cor_G2M = round(as.numeric(cor_g2m), 3)
)
cor_df
# 4. Harmony
combined <- RunHarmony(combined, group.by.vars = "orig.ident",
                       reduction = "pca", dims.use = 1:30)

# 5. 聚类 + UMAP
combined <- FindNeighbors(combined, dims = 1:30, reduction = "harmony")
combined <- FindClusters(combined, resolution = 0.3)
combined <- RunUMAP(combined, dims = 1:30, reduction = "harmony")