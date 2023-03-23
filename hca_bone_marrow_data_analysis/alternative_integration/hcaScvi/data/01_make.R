library(reticulate)
library(anndata)
library(sceasy)
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
library(Seurat)

ifnb <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/proc/integrate/ser/umap.rds')
ifnb <- ifnb@assays$RNA@counts
ifnb <- CreateSeuratObject(count=ifnb)

ifnb <- NormalizeData(ifnb, normalization.method = "LogNormalize", scale.factor = 10000)
ifnb <- FindVariableFeatures(ifnb, selection.method = "vst", nfeatures = 2000)
top2000 <- head(VariableFeatures(ifnb), 2000)
ifnb <- ifnb[top2000]
ifnb@meta.data$sample <- sub(':.*','',rownames(ifnb@meta.data))

adata <- convertFormat(ifnb, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata)

# run setup_anndata, use column stim for batch
scvi$data$setup_anndata(adata, batch_key = 'sample')

# create the model
model = scvi$model$SCVI(adata)

# train the model
model$train()

latent = model$get_latent_representation()

latent <- as.matrix(latent)
rownames(latent) = colnames(ifnb)
saveRDS(latent,file='/home-4/zji4@jhu.edu/scratch/diffpt/tb/latent.rds')



# ifnb[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(ifnb))
# 
# library(cowplot)
# # for jupyter notebook
# options(repr.plot.width=10, repr.plot.height=8)
# 
# ifnb <- RunUMAP(ifnb, dims = 1:10, reduction = "scvi", n.components = 2)
# p1 <- DimPlot(ifnb, reduction = "umap", group.by = "stim", pt.size=2)
# plot_grid(p1)
