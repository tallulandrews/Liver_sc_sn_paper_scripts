require("Seurat")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")


#if (!file.exists("All_merged_universal_genes_with_cluster.rds")) {
#	merged_obj <- readRDS("All_merged_universal_genes_harmony_integrated_v2.rds")

#	set.seed(3921)
#	merged_obj <- FindNeighbors(merged_obj, dims = 1:20, reduction="harmony")
#	set.seed(02789)
#	merged_obj <- FindClusters(merged_obj, resolution = 0.5)
#
#	saveRDS(merged_obj, "All_merged_universal_genes_with_cluster.rds")

#} else {
#	merged_obj <- readRDS("All_merged_universal_genes_with_cluster.rds")
#}

merged_obj <- readRDS("Sc_vs_sn_integrated_harmony_plus_analysis_Anno.rds")

#DimPlot(merged_obj, reduction="umap", group.by="sample", pt.size=0.1)

#DimPlot(merged_obj, reduction="umap", pt.size=0.1)

#Type_DimPlot(merged_obj, reduction="umap", type_col="marker_labs", cluster_col="marker_labs")

#anno_clusters <- table(merged_obj@meta.data$seurat_clusters, merged_obj@meta.data$consistent_labs)
#anno_clusters2 <- table(merged_obj@meta.data$seurat_clusters, merged_obj@meta.data$marker_labs)

#clusters <- rbind(
#	c("0", "Hepatocyte"),
#	c("1", "PortalHep"),
#	c("2", "Hepatocyte"),
#	c("3", "Hepatocyte"),
#	c("4", "cvLSECs"),
#	c("5", "Macrophage"),
#	c("6", "NKTcells"),
#	c("7", "Stellate"),
#	c("8", "Hepatocyte"),
#	c("9", "Cholangiocyte"),
#	c("10", "CentralHep"),
#	c("11", "interHep"),
#	c("12", "PortalLSECs"),
#	c("13", "Hepatocyte"),
#	c("14", "cvLSECs"),
#	c("15", "Bcell"),
#	c("16", "Hepatocyte"),
#	c("17", "LSECs"),
#	c("18", "Eryth")
#	)


#type <- clusters[merged_obj@meta.data$seurat_clusters,2]

types <- merged_obj@meta.data$All_Integrated_Manual
merged_obj@meta.data$Aggr_Type <- types


#### Cluster DE ####

type <- commandArgs(trailingOnly=TRUE)
#cluster <- as.character(cluster);
#for (cluster in sort(unique(merged_obj@meta.data$seurat_clusters))){
for (type in sort(unique(merged_obj@meta.data$Aggr_Type))){
	de_list <- list()
	de_score_mat <- c();
	print("cluster")
	print(type)
	col_names <- c();
	for (donor in sort(unique(merged_obj@meta.data$donor))) {
		print("donor")
		print(donor)
		if (sum(merged_obj@meta.data$Aggr_Type == as.character(type) & merged_obj@meta.data$donor == donor) == 0) {next;}
		subset <- merged_obj[,merged_obj@meta.data$Aggr_Type == as.character(type) &
						merged_obj@meta.data$donor == donor ]
		freqs <- table(factor(subset@meta.data$assay_type, levels=c("single_nuc", "single_cell")))
		if (min(freqs) < 5) {next;}
		print(dim(subset))
		print(table(subset@meta.data$assay_type))
		subset<- Seurat::NormalizeData(subset)
		de_this <- Seurat::FindMarkers(subset, slot="data", test.use="t", 
					group.by="assay_type", ident.1="single_cell", ident.2="single_nuc", 
					logfc.threshold=0)
		de_this$score <- de_this$avg_logFC
		de_this$score[de_this$p_val_adj > 0.05] <- 0
		de_this$rank <- rank(de_this$score);
		de_list[[paste(type, donor, sep="_")]] <- de_this

		out <- de_this$score
		names(out) <- rownames(de_this);
	
		if (length(de_score_mat) == 0) {
			de_score_mat <- matrix(out, ncol=1)
			rownames(de_score_mat) <- names(out)
		} else {
			all_genes <- sort(unique(c(names(out), rownames(de_score_mat))))
			out <- out[match(all_genes, names(out))]
			out[is.na(out)] <- 0
			de_score_mat <- de_score_mat[match(all_genes, rownames(de_score_mat)),]
			de_score_mat[is.na(de_score_mat)] <- 0
			de_score_mat <- cbind(de_score_mat, out);
			rownames(de_score_mat) <- all_genes;
		}
		col_names <- c(col_names, donor)
	}
	#colnames(de_score_mat) <- sort(unique(merged_obj@meta.data$donor))
	colnames(de_score_mat) <- col_names
	saveRDS(de_list, file=paste(type, "sn_vs_sc_de.rds", sep="_"))
	saveRDS(de_score_mat, file=paste(type, "sn_vs_sc_de_scoresonly.rds", sep="_"))

}



### EdgeR doesn't work because using individually scaled datasets.
#### Pseudobulk ####
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

pseudo_bulks <- get_pseudobulk(merged_obj@assays$RNA@data, merged_obj$Aggr_Type, merged_obj$sample)
pseudo_bulks <- t(t(pseudo_bulks)/colSums(pseudo_bulks)*10000)
metadata <- strsplit(colnames(pseudo_bulks), "_")
group <- sapply(metadata, function(x){x[[1]]})
indi <- sapply(metadata, function(x){x[[2]]})
assay <- sapply(metadata, function(x){if(length(x) == 3) {x[[3]]} else{"SC"}})
assay[grep("ST", assay)] <- "SN"
assay[grep("RESEQ", assay)] <- "SC"

pseudobulk_out <- list();

for (type in sort(unique(merged_obj@meta.data$Aggr_Type))){
	tmp <- pseudo_bulks[,group==type]
	sc_vs_sn <- assay[group==type]
	out <- t(apply(tmp, 1, function(x) {
		freqs <- table(factor(sc_vs_sn, levels=c("SN", "SC")))
		if (min(freqs) < 5) {return(c(0,0,1))}
		res <- t.test(x[sc_vs_sn=="SC"], x[sc_vs_sn=="SN"])
		return(c(res$estimate, res$p.value))
	}))
	colnames(out) <- c("SC", "SN", "pvalue")
	pseudobulk_out[[type]] <- out;
}

out <- t(apply(pseudo_bulks, 1, function(x) {
	res <- t.test(x[assay=="SC"], x[assay=="SN"])
	return(c(res$estimate, res$p.value))
}))
colnames(out) <- c("SC", "SN", "pvalue")
pseudobulk_out[["All"]] <- out

saveRDS(pseudobulk_out, file="SN_vs_SC_Merged_Harmony_Pseudobulk_DE.rds")

## Pathway Analysis
pseudobulk_out <- readRDS("SN_vs_SC_Merged_Harmony_Pseudobulk_DE.rds")

# UTRs
source("~/R-Scripts/Ensembl_Stuff.R")
utrs <-  read.table("/cluster/projects/macparland/TA/ExternalData/Ensembl/UTRs_ensembl_hsap_7_July2020.txt", header=TRUE)
utrs$symbol <- General_Map(utrs[,1], in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")

utrs_5 <- aggregate(utrs$utr_5_len, by=list(gene=utrs$symbol), max, na.rm=T)
utrs_5[!is.finite(utrs_5[,2]),2] <- 0
utrs_3 <- aggregate(utrs$utr_3_len, by=list(gene=utrs$symbol), max, na.rm=T)
utrs_3[!is.finite(utrs_3[,2]),2] <- 0

utr_bins <- c(-1, 1, 100, 500, 10000)
tmp <- cut(utrs_5[,2], breaks=utr_bins); levels(tmp) <- c("none", "short", "medium", "long")
utrs_5 <- cbind(utrs_5, tmp)
utr_bins <- c(-1, 1, 500, 1000, 100000)
tmp <- cut(utrs_3[,2], breaks=utr_bins); levels(tmp) <- c("none", "short", "medium", "long")
utrs_3 <- cbind(utrs_3, tmp)

# MiRNA
miRNA <- read.delim("/cluster/projects/macparland/TA/ExternalData/Ensembl/miRNA_binding_sites_Ensembl_biomart_7July2020.txt", sep="\t")
miRNA$symbol <- General_Map(miRNA[,5], in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")

miCount <- aggregate(miRNA[,3], by=list(miRNA$symbol), length)
miCount <- miCount[miCount[,1] != "",]
mi_bins <- c(-1, 10, 40, 1000)
tmp <- cut(miCount[,2], breaks=mi_bins); levels(tmp) <- c("few", "some", "many")
miCount <- cbind(miCount, tmp)
colnames(miCount) <- c("gene", "nmirna", "mirna_bin")

#splicing
transcripts <- readRDS("/cluster/projects/macparland/TA/ExternalData/Ensembl/gene_vs_transcript_ensembl.rds")
transcripts$symbol <- General_Map(transcripts[,1], in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
transcripts <- transcripts[transcripts$symbol != "",]
t_counts <- aggregate(transcripts[,2], by=list(transcripts$symbol), length)
tc_bins <- c(-1, 1, 7, 100000)
tmp <- cut(t_counts[,2], breaks=tc_bins); levels(tmp) <- c("unique", "few", "many")
t_counts <- cbind(t_counts, tmp)
colnames(t_counts) <- c("gene", "ntrans", "trans_bin")

#other gene stats
gene_info <- read.table("/cluster/projects/macparland/TA/ExternalData/Ensembl/Ensemble_gene_stats_various.tsv", header=T)
gene_info$gene_biotype <- as.character(gene_info$gene_biotype)
gene_info <- gene_info[gene_info$hgnc_symbol != "",]

gene_type <- gene_info[,c("hgnc_symbol", "gene_biotype", "chromosome_name")]
gene_type <- unique(gene_type)

gene_gc <- gene_info[,c("hgnc_symbol", "percentage_gene_gc_content")]
gene_gc <- unique(gene_gc)
gene_gc <-  aggregate(gene_gc[,2], by=list(gene=gene_gc[,1]), mean)
gc_bins <- c(0, 40, 60, 100)
binned <- cut(gene_gc[,2], breaks=gc_bins)
levels(binned) <- c("low", "med", "high")
gene_gc <- cbind(gene_gc, binned)



transcript_length <- gene_info[,c("hgnc_symbol", "transcript_length")]
tmp <- aggregate(transcript_length[,2], by=list(gene=transcript_length$hgnc_symbol), mean)
transcript_length <- tmp
tl_bins <- quantile(transcript_length[,2], c(0, 0.25, 0.5, 0.75, 1))
binned <-cut(transcript_length[,2], breaks=tl_bins)
levels(binned) <- c("Q1", "Q2", "Q3", "Q4")
transcript_length <- cbind(transcript_length, binned)


# Housekeeping genes
HK_genes <- read.table("/cluster/projects/macparland/TA/ExternalData/Eisenberg_2013_HK_genes.txt")

# Ribosomal
ribo_genes1 <- read.delim("/cluster/projects/macparland/TA/ExternalData/RPG_ribosomal_protein_gene_database_Human.tsv", sep="\t", header=F)
ribo_genes1 <- ribo_genes1[2:nrow(ribo_genes1),2]
mt_ribo <- ribo_genes1[grep("^MRP", ribo_genes1)]
ribo_genes2 <- c(as.character(gene_gc$gene[grepl("^RPS", gene_gc$gene)]), 
		as.character(gene_gc$gene[grepl("^RPL", gene_gc$gene)]))

# Exon / Intron lengths
gene_lengths <- read.table("/cluster/projects/macparland/TA/ExternalData/Ensembl/ensembl_gene_exon_length.tsv", sep="\t", header=T)
lengths <- data.frame(gene=gene_lengths[,4], exon=gene_lengths[,2], intron=gene_lengths[,5]-gene_lengths[,2], n_exon=gene_lengths[,3])
tl_bins <- quantile(lengths[,2], c(0, 0.25, 0.5, 0.75, 1))
binned <-cut(lengths[,2], breaks=tl_bins)
levels(binned) <- c("Q1", "Q2", "Q3", "Q4")
lengths$exon_split <- binned
tl_bins <- quantile(lengths[,3], c(0, 0.25, 0.5, 0.75, 1))
tl_bins[2] <- 1;
tl_bins[1] <- -1;
binned <-cut(lengths[,3], breaks=tl_bins)
levels(binned) <- c("Q1", "Q2", "Q3", "Q4")
lengths$intron_split <- binned
tl_bins <- quantile(lengths[,4], c(0, 0.25, 0.5, 0.75, 1))
tl_bins[1] <- 0;
binned <-cut(lengths[,4], breaks=tl_bins)
levels(binned) <- c("Q1", "Q2", "Q3", "Q4")
lengths$nexon_split <- binned

# Mitochondrial
mito_proteins <- read.delim("/cluster/projects/macparland/TA/ExternalData/Broad_Human_mitogenes_Human.MitoCarta3.0_MitoGenes.csv", sep=",", header=T)
mito_proteins <- as.character(mito_proteins$Symbol)
mito_proteins <- mito_proteins[mito_proteins != ""]

gene_type$gene_biotype[gene_type$hgnc_symbol %in% HK_genes[,1]] <- "house_keeping"
gene_type$gene_biotype[gene_type$hgnc_symbol %in% mito_proteins] <- "mitochondrial_protein"
gene_type$gene_biotype[gene_type$hgnc_symbol %in% ribo_genes2] <- "ribosomal_protein"
gene_type$gene_biotype[gene_type$hgnc_symbol %in% mt_ribo] <- "mt_ribosomal_protein"
gene_type$gene_biotype[gene_type$chromosome_name == "MT"] <- "mt_genome"
colnames(gene_type) <- c("gene", "biotype", "chromosome")
gene_type <- gene_type[!duplicated(gene_type[,1]),]

# Data things: gene_type, utrs_5, utrs_3, tanscript_length, gene_gc, t_counts, miCount, 

# analysis
global <- pseudobulk_out[["All"]]
#global <- global[!grepl("MT-", rownames(global)),]
#global <- global[!grepl("RPS", rownames(global)),]
#global <- global[!grepl("RPL", rownames(global)),]

#global_diff <- global[,"SC"]-global[,"SN"]
global_diff <- log2(global[,"SC"])-log2(global[,"SN"])

match_genes <- function(x, mat) {
	if (!is.null(dim(x))) {
		mat <- mat[match(rownames(x), mat$gene),];
	} else {
		mat <- mat[match(names(x), mat$gene),];
	}
	mat[is.na(mat)] <- 0;
	return(mat);
}

quantile_bins <- function(x) {
	bins <- quantile(x, c(0, 0.25, 0.5, 0.75, 1));
	bins[1] <- bins[1] -1;
	bins[5] <- bins[5] +1;
	binned <- cut(x, breaks=bins, right=TRUE);
	levels(binned) <- c("Q1", "Q2", "Q3", "Q4");
	return(binned)
}

global_diff <- global_diff[names(global_diff) %in% gene_type[,1]]

gene_type <- match_genes(global_diff, gene_type)
gene_type$biotype[grepl("IG_", gene_type$biotype)] <- "protein_coding"
gene_type$biotype[grepl("_pseudogene", gene_type$biotype)] <- "lncRNA"
gene_type$biotype <- factor(gene_type$biotype)

clean_diff <- global_diff[gene_type$biotype %in% c("house_keeping", "lncRNA", "protein_coding", "mitochondrial_protein")]

utrs_5_clean <- match_genes(clean_diff, utrs_5)
utrs_3_clean <- match_genes(clean_diff, utrs_3)
transcript_length_clean <- match_genes(clean_diff, transcript_length)
gene_gc_clean <- match_genes(clean_diff, gene_gc)
miCount_clean <- match_genes(clean_diff, miCount)
lengths_clean <- match_genes(clean_diff, lengths)

png("gene_stats_boxplots.png", width=6, height=8, units="in", res=300)
layout(rbind(c(1,1), c(2,3), c(4,5)))
# biotype
par(mar=c(6, 4, 1, 1))
gene_type$biotype <- factor(gene_type$biotype, levels=c(
		"protein_coding", 
		"house_keeping",
		"lncRNA",
		"ribosomal_protein",
		"mitochondrial_protein",
		"mt_ribosomal_protein",
		"mt_genome"))

boxplot(global_diff ~ gene_type$biotype, xlab="", ylab="SC-SN", notch=T, 
		names=c("protein coding", "housekeeping", "lncRNA", "ribosome protein", "mt protein", "mt ribosome", "mt genome"),
		col=c("#2c7fb8", "#41b6c4","#756bb1", "#7fcdbb", "#fed976", "#fd8d3c", "#f03b20"), 
		las=2)

par(mar=c(4,4,1,1))
# UTRs
utrs_5_clean[is.na(utrs_5_clean[,3]),3] <- "short"
utrs_5_clean[utrs_5_clean[,3]=="none",3] <- "short"
utrs_3_clean[is.na(utrs_3_clean[,3]),3] <- "short"
utrs_3_clean[utrs_3_clean[,3]=="none",3] <- "short"
boxplot(clean_diff ~ factor(utrs_5_clean[,3], levels = c("short", "medium", "long")),
	names=c("<100", "100-500", ">500"), 
	xlab="5' UTR (bp)", ylab="SC-SN",  notch=T, col=c("grey85", "grey65", "grey45"))
boxplot(clean_diff ~ factor(utrs_3_clean[,3], levels = c("short", "medium", "long")), 
	names=c("<500", "500-1000", ">1000"),
	xlab="3' UTR (bp)", ylab="SC-SN",  notch=T, col=c("grey85", "grey65", "grey45"))


#miRNA binding sites
boxplot(clean_diff ~ factor(miCount_clean[,3], levels=c("few", "some", "many")),
	names=c("< 10", "10-40", ">40"), 
	xlab="miRNA binding sites", ylab="SC-SN",  notch=T, col=c("grey85", "grey65", "grey45"))
# GC
boxplot(clean_diff ~ gene_gc_clean[,3], 
	names=c("< 40%", "40-60%", "> 60%"),
	xlab="GC content", ylab="SC-SN",  notch=T, 
	col=c("#998ec3", "grey85", "#f1a340"))
dev.off()


png("gene_stats_boxplots2.png", width=12*2/3, height=10, units="in", res=300)
par(mfrow=c(2,2))
par(mar=c(6,4,2,1))

#Transcript Length	
boxplot(clean_diff ~ quantile_bins(lengths_clean[,2]+lengths_clean[,3]), 
	names=c("Q1\n(< 15 kb)", "Q2\n(15 - 37 kb)", "Q3\n(37  - 86 kb)", "Q4\n(> 86 kb)"),
	main="Gene length", ylab="SC-SN",  notch=T, las=2, xlab="",
	col=c("grey85", "grey65", "grey45", "grey25"))

boxplot(clean_diff ~ quantile_bins(lengths_clean[,2]), 
	names=c("Q1\n(< 3.6 kb)", "Q2\n(3.6 - 5.5 kb)", "Q3\n(5.5 - 8 kb)", "Q4\n(> 8 kb)"),
	main="Exon length", ylab="SC-SN",  notch=T, las=2,xlab="",
	col=c("grey85", "grey65", "grey45", "grey25"))

boxplot(clean_diff ~ quantile_bins(lengths_clean[,3]), 
	names=c("Q1\n(< 10 kb)", "Q2\n(10kb - 30kb)", "Q3\n(30kb - 80kb)", "Q4\n(> 80kb)"),
	main="Inton length", ylab="SC-SN",  notch=T, las=2,xlab="",
	col=c("grey85", "grey65", "grey45", "grey25"))

boxplot(clean_diff ~ quantile_bins(lengths_clean[,4]), 
names=c("Q1\n(< 7)", "Q2\n(7 - 12)", "Q3\n(12 - 19)", "Q4\n(> 19)"),
	main="Number of Exons", ylab="SC-SN",  notch=T, las=2,xlab="",
	col=c("grey85", "grey65", "grey45", "grey25"))

dev.off()


#pch = biotype
keep <- rownames(global) %in% as.character(gene_lengths$symbol) & rownames(global) %in% as.character(gene_gc$gene)
global <- global[keep,]
gene_type <- gene_type[match(rownames(global), gene_type[,1]),]
pch = c(16, 16, 16, 17, 17, 16, 1)[gene_type$biotype]

gene_lengths <- gene_lengths[match(rownames(global), gene_lengths[,4]),]
gene_lengths$intron_length <- gene_lengths$genelength-gene_lengths$totalExonlength
bins = quantile(gene_lengths$intron_length, c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)/100)

#tmp2 <- log(gene_lengths$intron_length+1)
#tmp <- range(tmp2)
#bins <- seq(from=tmp[1]-1, to=tmp[2]+1, length=11)

bins[1] <- bins[1]-1;
cex = c(0.25, 0.5, 0.5, 0.8, 0.8, 1.2, 1.2, 2, 2, 3)[cut(gene_lengths$intron_length,bins)]
#cex = c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5)[cut(gene_lengths$intron_length,bins)]-0.5
#cex = c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5)[cut(tmp2,bins)]-1
#cex[cex <= 0] <- 0.5
#cex <- cex/2

gene_gc <- gene_gc[match(rownames(global), gene_gc[,1]),]
col=colorRampPalette(c("violet", "black"))(20)
tmp <- range(gene_gc[,2])
breaks <- seq(from=tmp[1]-1, to=tmp[2]+1, length=length(col)+1)
colours = col[cut(gene_gc[,2], breaks=breaks)]


png("gene_stats_scatter.png", width=8, height=8, units="in", res=300)
plot(log2(global[,"SC"]), log2(global[,"SN"]), pch=pch, cex=cex, col=colours, xlab="SC", ylab="SN")
dev.off()

#------------------------------------

gene_info <- gene_info[match(rownames(global), gene_info$hgnc_symbol),]
table(gene_info$gene_biotype)

mi_clean <- miCount[match(rownames(global), miCount[,1]),]
mi_clean[is.na(mi_clean[,2]),2] <- 0;

utrs_clean <- utrs[match(rownames(global), utrs$symbol),]
utrs_clean$utr_3_len[is.na(utrs_clean$utr_3_len)] <- 0
utrs_clean$utr_5_len[is.na(utrs_clean$utr_5_len)] <- 0

t_counts_clean <- t_counts[match(rownames(global), t_counts[,1]),]
t_counts_clean[is.na(t_counts[,2]),2] <- 1


png("diff_vs_stats.png", width=8, height=8, units="in", res=150)
par(mfrow=c(2,2))
plot(utrs_clean$utr_3_len, global_diff, pch=16, main="3' UTR", xlab="3' UTR length", ylab="SC-SN")
abline(h=0, col="red")
plot(utrs_clean$utr_5_len, global_diff, pch=16, main="5' UTR", xlab="5' UTR length", ylab="SC-SN")
abline(h=0, col="red")
plot(mi_clean[,2], global_diff, pch=16, main="miRNA binding sites", xlab="Num miRNA binding sites", ylab="SC-SN")
abline(h=0, col="red")
plot(t_counts_clean[,2], global_diff, pch=16, main="Num transcripts", xlab="Num transcripts", ylab="SC-SN")
abline(h=0, col="red")
dev.off()



