require("Seurat")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")

set.seed(3921)

prefix ="SN_SC"

# Which do we include in the integrated map?
dir <- "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned"
seurfiles <- c("C41_EmptyOnly.rds", 
	"C41_CST_EmptyOnly.rds",
	"C41_NST_EmptyOnly.rds",
	"C41_TST_EmptyOnly.rds",
	"C58_TST_EmptyOnly.rds",
	"C58_RESEQ_EmptyOnly.rds",
	"C70_TST_EmptyOnly.rds",
	"C70_RESEQ_EmptyOnly.rds",
	"C72_TST_EmptyOnly.rds",
	"C72_RESEQ_EmptyOnly.rds"
	);

samp_names <- unlist(lapply(strsplit(seurfiles, "_"), function(x){x <- x[c(-length(x))]; return(paste(x, collapse="_"))}))

obj_list <- list()
union_genes <- c();
for (i in 1:length(seurfiles)) {
	n <- samp_names[i];
	obj <- readRDS(paste(dir, seurfiles[i], sep="/"));

test <- c("CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "IGKC", "IGHE", "IGHM", "IGLC1", "IGLC2", "IGLC3")
print(sum(test %in% rownames(obj)))
	union_genes <- unique(c(union_genes, rownames(obj)))
	#Fix sample ID, and Donor ID
	obj@meta.data$sample <- obj@meta.data$orig.ident
	
	obj@meta.data$donor <- sapply(strsplit(as.character(obj@meta.data$sample), "_"), function(x){x[[1]]})

	# save sample specific clusters
	obj@meta.data$sample_specific_clusters <- paste(n, obj@meta.data$seurat_clusters, sep="_")

	# get rid of factors
	metadata_classes <- sapply(1:ncol(obj@meta.data), function(i){class(obj@meta.data[,i])})
	for (j in which(metadata_classes == "factor")) {
		obj@meta.data[,j] <- as.character(obj@meta.data[,j]);
	}

	
	obj <- Seurat::NormalizeData(obj, verbose = FALSE, normalization.method="LogNormalize", scale.factor=10000) 
	obj@meta.data$cell_barcode <- colnames(obj);
	obj@meta.data[[7]] <- obj@meta.data[[7]][[1]]
	obj@meta.data$sample <- rep(n, ncol(obj));
	obj@meta.data$cell_ID <- paste(obj@meta.data$sample, obj@meta.data$cell_barcode, sep="_")
	obj_list[[n]] <- obj
}

union_genes <- sort(union_genes)

hvgs <- c();
### union gene merged obj ###
require(Matrix)
all_Scaled <- c();
scaled_cell_ids <- c()
universal_genes <- c(-1)

for (i in 4:length(obj_list)) {
        n <- samp_names[i];
	obj <- obj_list[[i]];

	#obj_counts <- as.matrix(Seurat::GetAssayData(obj, "counts"))[match(union_genes, rownames(obj)), ]
	#obj_counts[is.na(obj_counts)] <- 0;
	#rownames(obj_counts) <- union_genes;

	#obj_counts <- Matrix::Matrix(obj_counts)
	#obj2 <- Seurat::CreateSeuratObject(counts=obj_counts, meta.data = obj@meta.data)
	obj <- Seurat::NormalizeData(obj, verbose = FALSE, normalization.method="LogNormalize", scale.factor=10000) 
	obj <- Seurat::ScaleData(obj, features=union_genes);
	scaled <- obj@assays$RNA@scale.data;
	scaled <- scaled[match(union_genes, rownames(scaled)),]
	scaled[is.na(scaled)] <- 0;
	scaled_cell_ids <- c(scaled_cell_ids, obj@meta.data$cell_ID);
	obj_list[[i]] <- obj;

	if (i == 1) {
		merged_obj <- obj
		all_Scaled <- scaled;
		universal_genes <- as.character(rownames(obj_list[[i]]))
		hvgs <- VariableFeatures(obj_list[[i]]);
		universal_genes <- as.character(rownames(obj_list[[i]]))
	} else {
		merged_obj2 <- merge(merged_obj, y=obj, add.cell.ids=c("", n), project="SC_SN_Map")
		scaled <- scaled[match(rownames(all_Scaled), rownames(scaled)),]
		all_Scaled <- cbind(all_Scaled, scaled);
		hvgs <- c(hvgs, VariableFeatures(obj_list[[i]]));
		universal_genes <- intersect(universal_genes, as.character(rownames(obj_list[[i]])))

dim(merged_obj);
dim(merged_obj2);
merged_obj <- merged_obj2;
	}
}
colnames(all_Scaled) <- scaled_cell_ids;
merged_obj@assays$RNA@scale.data <- all_Scaled

saveRDS(merged_obj, "SC_SN_AllGene_SeuratObject.rds")

# Keep HVGs seen in at least 2 datasets
hvgs <- unique(hvgs[duplicated(hvgs)])
hvgs <- hvgs[!grepl("^MT-", hvgs)]
hvgs <- hvgs[ hvgs %in% universal_genes]


# Merge Datasets
#### Merging does not merge individually scaled datasets!!

# Find common HVGs and detected genes
merged_obj <- NULL;
universal_genes <- c(-1)

hvgs <- c();
for (i in 1:length(obj_list)) {
        n <- samp_names[i];
	if (i == 1) {
		merged_obj <- obj_list[[i]]
		universal_genes <- as.character(rownames(obj_list[[i]]))
		hvgs <- VariableFeatures(obj_list[[i]]);
	} else {
		merged_obj <- merge(merged_obj, y=obj_list[[i]], add.cell.ids=c("", n), project="LiverMap")
		universal_genes <- intersect(universal_genes, as.character(rownames(obj_list[[i]])))
		hvgs <- c(hvgs, VariableFeatures(obj_list[[i]]));
	}
}

fix_names <- paste(merged_obj@meta.data$orig.ident, merged_obj@meta.data$cell_barcode, sep="_")
merged_obj <- RenameCells(merged_obj, new.names=fix_names)

# Keep HVGs seen in at least 2 datasets
hvgs <- unique(hvgs[duplicated(hvgs)])
hvgs <- hvgs[!grepl("^MT-", hvgs)]
hvgs <- hvgs[ hvgs %in% universal_genes]


# Scale within datasets
all_Scaled <- c();
scaled_cell_ids <- c()
for (i in 1:length(obj_list)) {
      n <- samp_names[i];
	obj <- obj_list[[i]]
	obj <- Seurat::ScaleData(obj, features=hvgs);
	
	scaled <- obj@assays$RNA@scale.data;
	scaled_cell_ids <- c(scaled_cell_ids, obj_list[[i]]@meta.data$cell_ID);
	if (i == 1) {
		all_Scaled <- scaled;
	} else {
		scaled <- scaled[match(rownames(all_Scaled), rownames(scaled)),]
		all_Scaled <- cbind(all_Scaled, scaled);
	}
}

colnames(all_Scaled) <- scaled_cell_ids;
merged_obj@assays$RNA@scale.data <- all_Scaled

merged_obj@misc$universal_genes <- universal_genes;
merged_obj@misc$repeated_hvgs <- hvgs;
merged_obj@misc$creation_date <- date();
VariableFeatures(merged_obj) <- hvgs;
merged_obj@meta.data$seurat_clusters <- paste(merged_obj@meta.data$orig.ident,
		as.character(merged_obj@meta.data$seurat_clusters), sep="_")

merged_obj@meta.data$assay_type <- rep("single_cell", nrow(merged_obj@meta.data))
nuclei <- grepl("ST", merged_obj@meta.data$sample)
merged_obj@meta.data$assay_type[nuclei] <- "single_nuc"

saveRDS(merged_obj, paste(prefix, "merged_obj.rds", sep="_"))

## BCR & TCR dotplot
require(ggplot2)
#cluster_ids <- readRDS("Sc_vs_sn_integrated_harmony_plus_analysis_Anno.rds")@meta.data
#cluster_ids <- cluster_ids[match(merged_obj@meta.data$cell_ID, rownames(cluster_ids)),]
#exclude <- is.na(rownames(cluster_ids));

cluster_ids <- readRDS("Integrated_with_Subannotations.rds")@meta.data
cluster_ids <- cluster_ids[match(merged_obj@meta.data$cell_ID, rownames(cluster_ids)),]


merged_obj@meta.data <- cluster_ids
cols <- c("#F8766D", "#00BFC4")

rownames(merged_obj@meta.data) <- colnames(merged_obj)

#tmp_obj <- merged_obj[,merged_obj@meta.data$All_Integrated_Manual %in% c("NKTcell", "Bcells")]
tmp_obj <- merged_obj[,grepl("Lymph", merged_obj@meta.data$sub_annotation)]
saveRDS(tmp_obj, "TCR_BCR_tmp_obj.rds")
png("ForPoster_BCR_TCR.png", width=9.5, height=13, unit="in", res=300)
#DotPlot(tmp_obj, features= c("CD3D", "CD3E", "CD8A", "CD8B", "CD4", "TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "IGKC", "IGHE", "IGHM", "IGLC1", "IGLC2", "IGLC3", "CD79A", "CD79B"), split.by="assay_type",cols=cols,dot.scale=30 group.by="All_Integrated_Manual") + coord_flip()
DotPlot(tmp_obj, features= c("NKG7", "GNLY", "GzMB", "GZMK", "CD3D", "CD3E", "CD8A", "CD8B", "CD4", "TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "IGKC", "IGHE", "IGHM", "IGLC1", "IGLC2", "IGLC3", "CD79A", "CD79B"), split.by="assay_type",cols=cols,dot.scale=20, group.by="sub_annotation") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#TCR_BCR <- c("CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "IGKC", "IGHE", "IGHM", "IGLC1", "IGLC2", "IGLC3")

## Back to actual script ##
make_umaps <- function(seur_obj, tag, reduc="pca") {
	set.seed(9428)

	seur_obj <- RunPCA(seur_obj, pc.genes = hvgs, 
			npcs = 20, verbose = FALSE)
	seur_obj <- RunUMAP(seur_obj, reduction=reduc, dims = 1:15, verbose = FALSE)

	png(paste(prefix, tag, "sample_umap.png", sep="_"), width=9, height =6, units="in", res=300)
	print(DimPlot(seur_obj, reduction="umap", group.by="sample", pt.size=0.1))
	dev.off();

	png(paste(prefix, tag, "assay_umap.png", sep="_"), width=9, height =6, units="in", res=300)
	print(DimPlot(seur_obj, reduction="umap", group.by="assay_type", pt.size=0.1))
	dev.off();

	png(paste(prefix, tag, "umap_mark_autoanno.png", sep="_"), width=12, height =6, units="in", res=300)
	print(Type_DimPlot(seur_obj, reduction="umap", type_col="marker_labs", cluster_col="marker_labs"))
	dev.off();

	png(paste(prefix, tag, "umap_scmap_autoanno.png", sep="_"), width=12, height =6, units="in", res=300)
	print(Type_DimPlot(seur_obj, reduction="umap", type_col="consistent_labs", cluster_col="marker_labs"))
	dev.off();
	return(seur_obj)
}


# raw individually scaled
merged_obj <- merged_obj[rownames(merged_obj) %in% universal_genes,]
merged_obj <- make_umaps(merged_obj, "indi_scaled")

# rescale across datasets
obj <- Seurat::ScaleData(merged_obj, features=hvgs);
obj <- make_umaps(obj, "rescaled")

# harmony individually scaled
require("harmony")
set.seed(10131)
merged_obj <- RunHarmony(merged_obj, c("sample"), plot_convergence = TRUE)
merged_obj <- make_umaps(merged_obj, "indi_scaled_harmony", reduc="harmony")
saveRDS(merged_obj, paste(prefix, "universal_genes_harmony.rds", sep="_"));

# harmony rescaled
require("harmony")
set.seed(10131)
obj <- RunHarmony(obj, c("sample"), plot_convergence = TRUE)
obj <- make_umaps(obj, "rescaled_harmony", reduc="harmony")


