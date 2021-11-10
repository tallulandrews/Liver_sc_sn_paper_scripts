obj <- readRDS("Sc_vs_sn_integrated_harmony_plus_analysis_Anno.rds")
require(Seurat)

type_barplot <- table(obj@meta.data$Manual_Core, obj@meta.data$sample)
type_barplot <- t(t(type_barplot)/colSums(type_barplot))

require(matrixStats)
sc_means <- rowMeans(type_barplot[,c(1,4,6,8)])
sc_vars <- rowVars(type_barplot[,c(1,4,6,8)])
sn_means <- rowMeans(type_barplot[,c(2,3,4,6,8,10)])
sn_vars <- rowVars(type_barplot[,c(2,3,4,6,8,10)])
sn_stderr <- sqrt(sn_vars)/sqrt(6)
sc_stderr <- sqrt(sc_vars)/sqrt(4)

cols <- c("#F8766D", "#00BFC4")

data <- cbind(sc_means, sn_means, sc_stderr, sn_stderr)
rownames(data) <- names(sc_means)

data <- data[order(data[,1], decreasing=T),]

pdf("Cell_type_Figure1.pdf", width=10, height=5)
barplot_obj <- barplot(t(data[,1:2])*100, beside=T, col=cols, 
				names=rownames(data), ylab="Frequency (%)", 
				ylim=c(0,100))

sc_top_errbar <- data[,1]+data[,3]*1.96
sn_top_errbar <- data[,2]+data[,4]*1.96

arrows(barplot_obj[1,], data[,1]*100, 
	barplot_obj[1,], sc_top_errbar*100, angle=90, len=0.1)

arrows(barplot_obj[2,], data[,2]*100, 
	barplot_obj[2,], sn_top_errbar*100, angle=90, len=0.1)

text(x=barplot_obj[1,], y=sc_top_errbar*100, round(data[,1]*100, digits=1), pos=3)
text(x=barplot_obj[2,], y=sn_top_errbar*100, round(data[,2]*100, digits=1), pos=3)

dev.off()