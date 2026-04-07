library(Seurat)
library(ggplot2)

#baysor output
seg <- read_csv("segmentation.csv")

#filter very low confidence assignments
tmp <- seg[seg$confidence > .5 & seg$assignment_confidence > .5 & !is.na(seg$cell), ]
mat <- table(tmp$gene, tmp$cell)
mat <- matrix(mat, ncol = ncol(mat), dimnames = dimnames(mat))

#baysor output
stats <- read_csv("segmentation_cell_stats.csv")
stats <- as.data.frame(stats)
rownames(stats) <- stats$cell

reg <- CreateSeuratObject(counts = mat, assay = "MERSCOPE", meta.data = as.data.frame(stats))

#remove low count cells 
reg <- reg[, reg$nCount_MERSCOPE > 10]

cents <- CreateCentroids(stats[Cells(reg), c("x", "y")])
cents@cells <- Cells(reg)
coords <- CreateFOV(coords =list(centroids = cents) ,
                    type = c("centroids"), 
                    molecules = NULL,
                    assay = "MERSCOPE")

reg[["REGION1"]] <- coords

reg <- SCTransform(reg, assay = "MERSCOPE", clip.range = c(-10, 10))
reg <- RunPCA(reg, npcs = 30, features = rownames(reg))
reg <- RunUMAP(reg, dims = 1:30)
reg <- FindNeighbors(reg, reduction = "pca", dims = 1:30)
reg <- FindClusters(reg, resolution = .7)

DimPlot(reg, label = T, repel=T, raster = F)
ImageDimPlot(reg)

FeaturePlot(reg, "RET", raster=F, label=T, repel=T)
FeaturePlot(reg, "S100B", raster=F, label=T, repel=T)
FeaturePlot(reg, "PTPRC", raster=F, label=T, repel=T)
FeaturePlot(reg, "CHGB", raster=F, label=T, repel=T)
FeaturePlot(reg, "EPCAM", raster=F, label=T, repel=T)
FeaturePlot(reg, "CD3D", raster=F, label=T, repel=T)
FeaturePlot(reg, "MS4A1", raster=F, label=T, repel=T)
FeaturePlot(reg, "PECAM1", raster=F, label=T, repel=T)
FeaturePlot(reg, "F3", raster=F, label=T, repel=T)
FeaturePlot(reg, "ADAMDEC1", raster=F, label=T, repel=T)
FeaturePlot(reg, "JCHAIN", raster=F, label=T, repel=T)
FeaturePlot(reg, "CCL19", raster=F, label=T, repel=T)
FeaturePlot(reg, "CCL21", raster=F, label=T, repel=T)
FeaturePlot(reg, "CD14", raster=F, label=T, repel=T)
FeaturePlot(reg, "PROX1", raster=F, label=T, repel=T)
FeaturePlot(reg, "KIT", raster=F, label=T, repel=T)
FeaturePlot(reg, "ANO1", raster=F, label=T, repel=T)
FeaturePlot(reg, "C3", raster=F, label=T, repel=T)
FeaturePlot(reg, "C7", raster=F, label=T, repel=T)

reg <- RenameIdents(reg, "13"="Stromal 1", "12"="Neurons&Glia",  "0"="Epithelium",
                     "3"="Epithelium", "4"="Epithelium", "7"="T-cells", "21"="Epithelium",
                     "14"="Epithelium","15"="Epithelium", "23"="Epithelium", "9"="Stromal 4",
                     "11"="Epithelium","2"="Stromal 3",
                     "20"="ICCs", "19"="Lymphatics", "16"="Myeloid", "6"="B-cells",
                     "5"="B-cells","10"="Endothelium",  "17"="Stromal 2", "18"="T-cells",
                     "8"="Stromal 3", "22"="Mast",
                     "1"="Muscularis")

saveRDS(reg, file="seurat.RDS")

