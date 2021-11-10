require(Seurat)
source("~/scripts/LiverMap2.0/My_R_Scripts.R") # pseudobulk
source("~/scripts/LiverMap2.0/Setup_autoannotation.R") # matching/labelling clusters

# Read in matched files

obj_sn = readRDS("C41_TST_EmptyOnly.rds")
obj_sc = readRDS("C41_EmptyOnly.rds")

# Match clusters - no integration b/c differences between clusters more conserved than overall similarity.
# ---- This is based on scaled average expression relative to other clusters or DE genes
# ---- One option: enrichment for dataset A's markers in dataset B's markers
sn_scmap_cluster_anno <- cell_anno_to_cluster_anno(obj_sn@meta.data$consistent_labs, obj_sn@meta.data$seurat_clusters)
sn_cluster_anno <- cell_anno_to_cluster_anno(obj_sn@meta.data$marker_labs, obj_sn@meta.data$seurat_clusters)
sn_cluster_labs <- as.character(sn_scmap_cluster_anno[,2])
replace <- grepl("Hep", sn_cluster_labs) | sn_cluster_labs == "ambiguous"
sn_cluster_labs[replace] <- as.character(sn_cluster_anno[replace,2])


sc_scmap_cluster_anno <- cell_anno_to_cluster_anno(obj_sc@meta.data$consistent_labs, obj_sc@meta.data$seurat_clusters)
sc_cluster_anno <- cell_anno_to_cluster_anno(obj_sc@meta.data$marker_labs, obj_sc@meta.data$seurat_clusters)
sc_cluster_labs <- as.character(sc_scmap_cluster_anno[,2])
replace <- grepl("Hep", sc_cluster_labs) | sc_cluster_labs == "ambiguous"
sc_cluster_labs[replace] <- as.character(sc_cluster_anno[replace,2])


hvgs <- intersect(VariableFeatures(obj_sn), VariableFeatures(obj_sc))
hvgs <- hvgs[hvgs %in% rownames(obj_sn) & hvgs %in% rownames(obj_sc)]
sc_pseudo <- group_rowmeans(obj_sc@assays$RNA@data, obj_sc@meta.data$seurat_clusters, type="sum") 
sc_pseudo <- t(t(sc_pseudo)/Matrix::colSums(sc_pseudo)*10000)
sn_pseudo <- group_rowmeans(obj_sn@assays$RNA@data, obj_sn@meta.data$seurat_clusters, type="sum") 
sn_pseudo <- t(t(sn_pseudo)/Matrix::colSums(sn_pseudo)*10000)
sc_pseudo <- sc_pseudo[rownames(sc_pseudo) %in% hvgs,]
sn_pseudo <- sn_pseudo[rownames(sn_pseudo) %in% hvgs,]

require(proxy)
similar<-proxy::simil(t(sc_pseudo), t(sn_pseudo), method="cosine")
tmp <- group_rowmeans(similar, colnames(similar), type="mean")
tmp <- group_colmeans(tmp, rownames(similar), type="mean")
colnames(similar) <- sn_cluster_labs
rownames(similar) <- sc_cluster_labs

test <- t(apply(similar, 1, function(x){x[x < max(x)]<-0; return(x)}))
colnames(test) <- sn_cluster_labs
rownames(test) <- sc_cluster_labs

test2 <- apply(similar, 2, function(x){x[x < max(x)]<-0; return(x)})
colnames(test2) <- sn_cluster_labs
rownames(test2) <- sc_cluster_labs
test*test2

# I've changed my mind, I think integrating all of them together is the best way to go, 
# however may want to modify scaling to make this easier. <- need to match clusters across 
# all of the datasets, not just the pairs to look at consistency of markers trends 
# across individuals

# Within each cluster 
# -- DE between sn and sc
# wilcox.test?
# t.test?
for (c in matched_clusters) {
	

}


# Find consistent across samples -> use script already working on to clean up patient DE from background.

# Examine consistent differences for trends w.r.t. gene length, splicing, miRNA, UTR, etc..

# Also look at reciprocality of marker enrichments. Are sn markers more enriched in sc markers or sc markers more enriched in sn markers?
# ------- does this make sense? Are these symetric?
# ------- I think I need to use GSEA for this...

source("/cluster/home/tandrews/R-Scripts/GSEA_code.R")
