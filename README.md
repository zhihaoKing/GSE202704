scRNA-Seq microglia gse202704
==========
## Data preparing
```r
library(data.table)
library(Matrix)

s <- c("C17H","C17P","C18H","C18P","C19H","C19P","G5H","G5P","G6H","G6P","G7H","G7P",
       "S30H","S30P","S31H","S31P","S32H","S32P")
m <- list()
b <- list()
g <- list()
matrix <- list()
barcodes <- list()
gene <- list()

for (i in seq(s)) {
  m[[i]] <- paste0(s[i],".matrix.mtx")
  b[[i]] <- paste0(s[i],".barcodes.tsv")
  g[[i]] <- paste0(s[i],".genes.tsv")
}

for (i in seq(s)) {
  matrix[[i]] <- readMM(file = m[[i]])
  barcodes[[i]] <- fread(b[[i]],header = F)
  gene[[i]] <- fread(g[[i]],header = F)
}

for (i in seq(s)) {
  colnames(matrix[[i]]) <- barcodes[[i]]$V1
  rownames(matrix[[i]]) <- gene[[i]]$V2
}

for (i in seq(s)) {
  barcodes[[i]]$sample <- s[i]
}

barcodes[[1]]$group <- "CGF"
barcodes[[2]]$group <- "CGF"
barcodes[[3]]$group <- "CGF"
barcodes[[4]]$group <- "CGF"
barcodes[[5]]$group <- "CGF"
barcodes[[6]]$group <- "CGF"
barcodes[[7]]$group <- "GF"
barcodes[[8]]$group <- "GF"
barcodes[[9]]$group <- "GF"
barcodes[[10]]$group <- "GF"
barcodes[[11]]$group <- "GF"
barcodes[[12]]$group <- "GF"
barcodes[[13]]$group <- "SPF"
barcodes[[14]]$group <- "SPF"
barcodes[[15]]$group <- "SPF"
barcodes[[16]]$group <- "SPF"
barcodes[[17]]$group <- "SPF"
barcodes[[18]]$group <- "SPF"
meta <- rbind(barcodes[[1]],barcodes[[2]],barcodes[[3]],barcodes[[4]],
              barcodes[[5]],barcodes[[6]],barcodes[[7]],barcodes[[8]],
              barcodes[[9]],barcodes[[10]],barcodes[[11]],barcodes[[12]],
              barcodes[[13]],barcodes[[14]],barcodes[[15]],barcodes[[16]],
              barcodes[[17]],barcodes[[18]])
#setDT(meta)[ , ID := .GRP, by = rownames(meta)] 
ex <- cbind(matrix[[1]],matrix[[2]],matrix[[3]],matrix[[4]],matrix[[5]],
            matrix[[6]],matrix[[7]],matrix[[8]],matrix[[9]],matrix[[10]],
            matrix[[11]],matrix[[12]],matrix[[13]],matrix[[14]],matrix[[15]],
            matrix[[16]],matrix[[17]],matrix[[18]])
colnames(ex) <- rownames(meta)
save(meta,ex,file="gse202704.RData")
```
## Seurat Object
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(Seurat)
library(ggplot2)
library(tidyr)

setwd("/home/shpc_101115/GSE202704")
load("/home/shpc_101115/GSE202704/gse202704.RData")#load expression profile and meta information

glia <- CreateSeuratObject(counts = ex,
                            meta.data = meta,
                            min.cells = 3, 
                            min.features = 50)##creat Seurat Object

glia@meta.data$tissue <- ifelse(glia@meta.data$sample=="C17H"|glia@meta.data$sample=="C18H"|
                                glia@meta.data$sample=="C19H"|glia@meta.data$sample=="G5H"|
                                glia@meta.data$sample=="G6H"|glia@meta.data$sample=="G7H"|
                                glia@meta.data$sample=="S30H"|glia@meta.data$sample=="S31H"|
                                glia@meta.data$sample=="S32H","Hippocampus","Prefrontal cortex")

#save(glia,file="glia_seuratobject.RData")

#load("~/GSE202704/glia_seuratobject.RData")
cells.use <- row.names(glia@meta.data)[which(glia@meta.data$tissue=='Hippocampus')]
glia_h <-subset(glia, cells=cells.use)  
save(glia_h,file = "/home/shpc_101115/GSE202704/Hip/glia_h.RData")

cells.use <- row.names(glia@meta.data)[which(glia@meta.data$tissue=='Prefrontal cortex')]
glia_p <-subset(glia, cells=cells.use)  
save(glia_p,file = "/home/shpc_101115/GSE202704/PFC/glia_p.RData")
```
## Analyze glia_p 
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(Seurat)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel) 

load("~/GSE202704/Hip/glia_h.RData")
```
## QC and Normalize data and feature selection and Scaling the data and Perform linear/non-linear dimensional reduction
```r
glia_h[["percent.mt"]] <- PercentageFeatureSet(glia_h, pattern = "^MT-")

glia_h <- subset(glia_h, subset = nFeature_RNA > 200 & percent.mt < 10)
glia_h <- NormalizeData(glia_h, normalization.method = "LogNormalize", scale.factor = 10000)
glia_h <- FindVariableFeatures(glia_h, selection.method = "vst", nfeatures = 2000)

#head(VariableFeatures(glia_h), 10)

all.genes <- rownames(glia_h)
glia_h <- ScaleData(glia_h, features = all.genes)
glia_h <- RunPCA(glia_h, features = VariableFeatures(object = glia_h))

glia_h <- FindNeighbors(glia_h, dims = 1:10)
glia_h <- FindClusters(glia_h, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2))

Idents(glia_h) <- glia_h@meta.data$RNA_snn_res.0.8
glia_h <- RunUMAP(glia_h, dims = 1:10)

pdf("/home/shpc_101115/GSE202704/Hip/glia_h_umap.pdf",w=10,h=10)
DimPlot(glia_h, reduction = "umap",label = T,label.size = 3.6)
dev.off()
```
## Cell annotation and select microglia
```r
genes_to_check = c("Grin2a","Syt1","Grin1","Gad1","Gad2","Plp1","Mog","Mbp",
                   "Pdgfra","Vcan","Csf1r","Ctss","C1qa","Aqp4","Gja1","Pecam1",
                   "Cldn5","Kdr","Vwf","Dcn","Col1a1")

p <- DotPlot(glia_h, features = genes_to_check) + coord_flip()

pdf("/home/shpc_101115/GSE202704/Hip/glia_h_check_gene.pdf",w=10,h=10)
p
dev.off()

save(glia_h,file = "/home/shpc_101115/GSE202704/Hip/glia_h1.RData")
load("/home/shpc_101115/GSE202704/Hip/glia_h1.RData")
#select microgila

microglia_h <- glia_h[,glia_h@meta.data$RNA_snn_res.0.8 %in% c(9,34)]
microglia_h@meta.data <- microglia_h@meta.data[,c(1:8)]
Idents(microglia_h) <- microglia_h@meta.data$group
#microglia_h <- RunPCA(microglia_h, features = rownames(microglia_h))
#microglia_h <- RunUMAP(microglia_h, features = rownames(microglia_h))
#pdf("/home/shpc_101115/GSE202704/Hip/pca_microglia_h.pdf",w=4,h=4)
DimPlot(microglia_h, reduction = "pca",dims = c(1,2))
#dev.off()

save(microglia_h,file = "/home/shpc_101115/GSE202704/Hip/microglia_h.RData")
```
## Difference Analysis
```r
h_de_gf_cgf <- FindMarkers(microglia_h, ident.1 = "GF", ident.2 = "CGF",test.use = "wilcox")

### 差异分析火山图 ###

h_de_gf_cgf$gene <- rownames(h_de_gf_cgf)

h_de_gf_cgf <- h_de_gf_cgf %>% 
  mutate(Expression = case_when(avg_log2FC >= log(2) & p_val_adj <= 0.05 ~ "Up-regulated",
                           avg_log2FC <= -log(2) & p_val_adj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged"))

p4 <- ggplot(h_de_gf_cgf, aes(avg_log2FC, -log(p_val_adj,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

top <- 10
top_genes <- bind_rows(
  h_de_gf_cgf %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>% 
    head(top),
  h_de_gf_cgf %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>% 
    head(top)
)

p5 <-  p4 +
  geom_label_repel(data = top_genes,
                   mapping = aes(avg_log2FC, -log(p_val_adj,10), label = gene),
                   size = 1.8)

pdf("/home/shpc_101115/GSE202704/Hip/microglia_de_GF_CGF.pdf",w=5,h=3.5)
p5
dev.off()
save(h_de_gf_spf,h_de_gf_cgf,file = "/home/shpc_101115/GSE202704/Hip/Hip_DE_result.Rdata")
```
