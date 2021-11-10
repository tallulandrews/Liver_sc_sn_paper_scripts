#matched_files <- list(
#	Tcell=c("SC_Empty_harmony_CD3abTcells_DE.rds", "SN_Empty_harmony_Tcell_DE.rds"), 
#	NKcell=c("SC_Empty_harmony_NKcells_DE.rds", "SN_Empty_harmony_Tcell_DE.rds"), 
#	Stellate=c("SC_Empty_harmony_Stellate_DE.rds", "SN_Empty_harmony_Stellate_DE.rds"), 
#	PortHep=c("SC_Empty_harmony_PortalHep_DE.rds", "SN_Empty_harmony_PortalHep_DE.rds"), 
#	InterHep=c("SC_Empty_harmony_InterHep_DE.rds", "SN_Empty_harmony_InterHep_DE.rds"), 
#	PortEndo=c("SC_Empty_harmony_PortalEndo_DE.rds", "SN_Empty_harmony_PortalEndo_DE.rds"),
#	cvLSECs=c("SC_Empty_harmony_cvLSECs_DE.rds", "SN_Empty_harmony_cvLSECs_DE.rds"),
#	Cholan=c("SC_Empty_harmony_Cholangiocyte_DE.rds", "SN_Empty_harmony_Cholangiocyte_DE.rds"),
#	CentHep=c("SC_Empty_harmony_CentralHep_DE.rds", "SN_Empty_harmony_CentralHep_DE.rds"),
#	InfMac=c("SC_Empty_harmony_InfMac_DE.rds", "SN_Empty_harmony_NonInfMac_DE.rds"),
#	NonInfMac1=c("SC_Empty_harmony_NonInfMac1_DE.rds", "SN_Empty_harmony_NonInfMac_DE.rds"),
#	NonInfMac2=c("SC_Empty_harmony_NonInfMac2_DE.rds", "SN_Empty_harmony_NonInfMac_DE.rds"))

matched_files <- list(
	NKTcell=c("SC_All_Manual_Core_Celltype_Markers_NKTcell_DE.rds", 
		"SN_All_Manual_Core_Celltype_Markers_NKTcell_DE.rds"), 
	Macro=c("SC_All_Manual_Core_Celltype_Markers_Macrophage_DE.rds", 
		"SN_All_Manual_Core_Celltype_Markers_Macrophage_DE.rds"), 
	Bcell=c("SC_All_Manual_Core_Celltype_Markers_Bcells_DE.rds", 
		"SN_All_Manual_Core_Celltype_Markers_Bcells_DE.rds"), 
	Cholan=c("SC_All_Manual_Core_Celltype_Markers_Cholangiocyte_DE.rds", 
		"SN_All_Manual_Core_Celltype_Markers_Cholangiocyte_DE.rds"), 
	Hepato=c("SC_All_Manual_Core_Celltype_Markers_Hepatocyte_DE.rds", 
		"SN_All_Manual_Core_Celltype_Markers_Hepatocyte_DE.rds"), 
	InfMac=c("SC_All_Manual_Core_Celltype_Markers_InfMac_DE.rds", 
		"SN_All_Manual_Core_Celltype_Markers_InfMac_DE.rds"), 
	NonInfMac=c("SC_All_Manual_Core_Celltype_Markers_NonInfMac_DE.rds", 
		"SN_All_Manual_Core_Celltype_Markers_NonInfMac_DE.rds"), 
	PortalEndo=c("SC_All_Manual_Core_Celltype_Markers_PortalEndo_DE.rds", 
		"SN_All_Manual_Core_Celltype_Markers_PortalEndo_DE.rds"), 
	Stellate=c("SC_All_Manual_Core_Celltype_Markers_Stellate_DE.rds", 
		"SN_All_Manual_Core_Celltype_Markers_Stellate_DE.rds"), 
	LSECs=c("SC_All_Manual_Core_Celltype_Markers_LSECs_DE.rds", 
		"SN_All_Manual_Core_Celltype_Markers_LSECs_DE.rds"))
tag = "All_Manual_Harmony"

for(id in names(matched_files)) {
scDE <- readRDS(matched_files[[id]][1])
snDE <- readRDS(matched_files[[id]][2])


scDE <- scDE$seur_wilcox
snDE <- snDE$seur_wilcox


common_genes <- sort(intersect(rownames(scDE), rownames(snDE)))

common_genes <- common_genes[!grepl("MT-", common_genes)]

png(paste(tag, id, "SC_vs_SN_DE.png", sep="_"), width=8, height=5, units="in", res=300)
layout(mat=rbind(c(1,3), c(2,3)), widths=c(1,2), heights=c(1,1))
require(venn)
venn(x=list(rownames(scDE)[scDE$avg_logFC > 0 & scDE$p_val_adj < 0.05], 
	    rownames(snDE)[snDE$avg_logFC > 0 & snDE$p_val_adj < 0.05]),
		snames=c("Cells", "Nuclei"), box=F, zcolor=c("red", "yellow"))
title(main="Positive Markers", line=-2.5)

venn(x=list(rownames(scDE)[scDE$avg_logFC < 0 & scDE$p_val_adj < 0.05], 
	    rownames(snDE)[snDE$avg_logFC < 0 & snDE$p_val_adj < 0.05]),
		snames=c("Cells", "Nuclei"), box=F, zcolor=c("blue", "violet"))
title(main="Negative Markers", line=-2.5)

scDE_reorder <- scDE[match(common_genes, rownames(scDE)),]
snDE_reorder <- snDE[match(common_genes, rownames(snDE)),]

#reg=lm(snDE_reorder$avg_logFC~scDE_reorder$avg_logFC)
#tmp <- summary(reg)
reg <- list(coefficients=c(0,1), residuals=snDE_reorder$avg_logFC-scDE_reorder$avg_logFC)

best <- snDE_reorder$avg_logFC+scDE_reorder$avg_logFC
# collect genes to label
key_genes <- reg$residuals <= quantile(reg$residuals, prob=10/length(common_genes)) |
	reg$residuals >= quantile(reg$residuals, prob=1-10/length(common_genes)) |
	best <= quantile(best, prob=5/length(common_genes)) |
	best >= quantile(best, prob=1-20/length(common_genes))

#key_genes <- c(common_genes[abs(reg$residuals) > 2.5], common_genes[abs(best) > 2.5] )
key_genes <- common_genes[key_genes]
labls <- common_genes;
labls[!labls %in% key_genes] <- ""
labl_pos <- rep(0, length(common_genes));  # 1=down, 2=left, 3=up, 4=right
labl_pos[scDE_reorder$avg_logFC > 0 & reg$residuals >0] <- 2
labl_pos[scDE_reorder$avg_logFC > 0 & reg$residuals <0] <- 1
labl_pos[scDE_reorder$avg_logFC < 0 & reg$residuals >0] <- 3
labl_pos[scDE_reorder$avg_logFC < 0 & reg$residuals <0] <- 4

#Colour by signficance!
dot_col <- rep("grey75", length(common_genes));
dot_col[scDE_reorder$p_val_adj < 0.05 & snDE_reorder$p_val_adj < 0.05] <- "grey35"
dot_col[common_genes %in% key_genes] <- "black"

par(mar=c(4,4,1,1), xpd=FALSE)
plot(scDE_reorder$avg_logFC, snDE_reorder$avg_logFC, xlab="Single Cell (log2FC)", ylab="Single Nuc (log2FFC)", pch=16, col=dot_col, cex=0.5)

text(scDE_reorder$avg_logFC, snDE_reorder$avg_logFC, labels=labls, cex=0.5, pos=labl_pos)

abline(a = reg$coefficients[1], b=reg$coefficients[2], col="red", lwd=2)
#paste("R2 =",signif(tmp$adj.r.squared, digits=1))
legend("topleft", c("SN=SC", paste("Spearman =", 
		signif(cor(scDE_reorder$avg_logFC,snDE_reorder$avg_logFC, method="spearman"), 
		digits=2))), lwd=2, lty=1, 
		col=c("red", "white"), bty="n")

dev.off()

out_tab <- data.frame(Gene=common_genes, sc_l2fc=scDE_reorder$avg_logFC, sc_qval=scDE_reorder$p_val_adj, sn_l2fc=snDE_reorder$avg_logFC, sn_qval=snDE_reorder$p_val_adj)
write.table(out_tab, file=paste(tag, id, "SC_vs_SN_gene_list.csv", sep="_"), sep=",")

}
