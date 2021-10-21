library(Seurat)
library('biomaRt')
library(fgsea)
library(ggplot2)
library(gplots)
library(data.table)
library(RColorBrewer)

library(dplyr)
library(Matrix)
library(cowplot)

bothdata<-readRDS("bothdata.rds")
Fibroblast <- t(as.matrix(GetAssayData(bothdata, slot = "count")[, WhichCells(bothdata, ident = c("Fibroblast_TUMOR","Fibroblast_CTRL"))]))

tmp1=0
matrix<-matrix(, nrow=nrow(Fibroblast),ncol=ncol(Fibroblast))
rownames(matrix)<-rownames(Fibroblast)
colnames(matrix)<-colnames(Fibroblast)
for (c in 1:ncol(Fibroblast)){
  tmp<-0
  for (r in 1:nrow(Fibroblast)){
    if (Fibroblast[r,c] >=1){
      tmp <- tmp+1
    }
    }
  if (tmp/nrow(Fibroblast) >=0.001){
  tmp1<-tmp1+1
  matrix[,tmp1]<-Fibroblast[,c]
  }
}
dat<-matrix[ , colSums(is.na(matrix)) == 0]


tf<-as.data.frame(dat)
cov<-sapply(tf, function(x) sd(x) / mean(x) * 100)
mcov<-as.matrix(cov)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
up100<-read.table("up100.txt",head=F)
up_gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values = up100, mart= mart)
up100gs=up_gene_IDs$hgnc_symbol

dn100<-read.table("dn100.txt",head=F)
dn_gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values = dn100, mart= mart)
dn100gs=dn_gene_IDs$hgnc_symbol


list_up100gs<-list(up100gs)
names(list_up100gs)<-c("up100")
a<-as.vector(mcov)
names(a)<-rownames(mcov)
mcov_sort = sort(a, decreasing = TRUE)
mcov_sort = mcov_sort[!duplicated(names(mcov_sort))]
mcov_sort<-mcov_sort[!is.na(mcov_sort) & !is.infinite(mcov_sort)]
fgsea_up <- fgsea(list_up100gs,mcov_sort,nperm=10000)
plotEnrichment(list_up100gs[["up100"]],mcov_sort) + labs(title="up 100")

list_dn100gs<-list(dn100gs)
names(list_dn100gs)<-c("dn100")
fgsea_dn <- fgsea(list_dn100gs,mcov_sort,nperm=10000)
plotEnrichment(list_dn100gs[["dn100"]],mcov_sort) + labs(title="dn 100")

##########
common_up <- intersect(up100gs, colnames(tf))
tf_up_100<-subset(tf,select=common_up)
cov_up_100<-sapply(tf_up_100, function(x) sd(x) / mean(x) * 100)
common_dn <- intersect(dn100gs, colnames(tf))
tf_dn_100<-subset(tf,select=common_dn)
cov_dn_100<-sapply(tf_dn_100, function(x) sd(x) / mean(x) * 100)
hist(cov_up_100,breaks=20)
hist(cov_dn_100,breaks=20)
dev.off()
########
Fibroblast <- subset(bothdata, idents = c("Fibroblast_TUMOR","Fibroblast_CTRL"))
Idents(object = Fibroblast) <- Fibroblast@meta.data$orig.ident
xlist<-paste("T",seq(1,24,by=1),sep="")
for (tumor in xlist) {
tryCatch(
{
Fibroblast.DE <- FindMarkers(Fibroblast, ident.1 = tumor, ident.2 = c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11"), verbose = FALSE, logfc.threshold=0, min.pct=0, min.cells.feature =1, min.cells.group =1)
write.table(Fibroblast.DE,paste(tumor,"_Fibroblast_DE.txt",sep=""),sep="\t")
},
error = function(e) {
}
)
}
########correlation
Fibroblast <- subset(bothdata, idents = c("Fibroblast_TUMOR","Fibroblast_CTRL"))
Idents(object = Fibroblast) <- Fibroblast@meta.data$orig.ident

matrix<-Fibroblast$RNA@data

matrix_mod<-as.matrix(matrix)
allgenes<-rownames(matrix_mod)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
up100<-read.table("up100.txt",head=F)
up_gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values = up100, mart= mart)
up100gs=up_gene_IDs$hgnc_symbol

dn100<-read.table("dn100.txt",head=F)
dn_gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values = dn100, mart= mart)
dn100gs=dn_gene_IDs$hgnc_symbol


common_up <- intersect(up100gs, allgenes)
common_dn <- intersect(dn100gs, allgenes)
cor_matrix_mod<-cor(t(matrix_mod))
gene_up<-matrix_mod[common_up,]
cor_up_matrix <- cor(t(gene_up))
gene_dn<-matrix_mod[common_dn,]
cor_dn_matrix <- cor(t(gene_dn))
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap.2(cor_up_matrix, dendrogram='none',main ='up_100', Rowv=F, Colv=F, trace='none', col = colorRampPalette(c('blue', 'yellow'))(12), density.info='none',keysize=1)
heatmap.2(cor_dn_matrix, dendrogram='none',main ='dn_100', Rowv=F, Colv=F, trace='none', col = colorRampPalette(c('blue', 'yellow'))(12), density.info='none',keysize=1)

heatmap.2(cor_matrix_mod, dendrogram='none', main ='all genes', Rowv=Null, Colv=Null, trace='none', col = colorRampPalette(c('blue', 'yellow'))(12), labRow=NA, labCol=NA , density.info='none',keysize=1)
####FILTER
bothdata<-readRDS("bothdata.rds")
Fibroblast <- t(as.matrix(GetAssayData(bothdata, slot = "data")[, WhichCells(bothdata, ident = c("Fibroblast_TUMOR","Fibroblast_CTRL"))]))

tmp1=0
matrix<-matrix(, nrow=nrow(Fibroblast),ncol=ncol(Fibroblast))
rownames(matrix)<-rownames(Fibroblast)
colnames(matrix)<-colnames(Fibroblast)
for (c in 1:ncol(Fibroblast)){
  tmp<-0
  for (r in 1:nrow(Fibroblast)){
    if (Fibroblast[r,c] >0.1){
      tmp <- tmp+1
    }
    }
  if (tmp/nrow(Fibroblast) >=0.05){
  tmp1<-tmp1+1
  matrix[,tmp1]<-Fibroblast[,c]
  }
}
dat<-matrix[ , colSums(is.na(matrix)) == 0]
matrix_mod<-as.matrix(t(dat))
allgenes<-rownames(matrix_mod)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
up100<-read.table("up100.txt",head=F)
up_gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values = up100, mart= mart)
up100gs=up_gene_IDs$hgnc_symbol

dn100<-read.table("dn100.txt",head=F)
dn_gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values = dn100, mart= mart)
dn100gs=dn_gene_IDs$hgnc_symbol


common_up <- intersect(up100gs, allgenes)
common_dn <- intersect(dn100gs, allgenes)

#cor_matrix_mod<-cor(t(matrix_mod))
gene_up<-matrix_mod[common_up,]
cor_up_matrix <- cor(t(gene_up))
gene_dn<-matrix_mod[common_dn,]
cor_dn_matrix <- cor(t(gene_dn))
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap.2(cor_up_matrix, dendrogram='none',main ='up_100', Rowv=F, Colv=F, trace='none', col = colorRampPalette(c('blue', 'yellow'))(12), density.info='none',keysize=1)
heatmap.2(cor_dn_matrix, dendrogram='none',main ='dn_100', Rowv=F, Colv=F, trace='none', col = colorRampPalette(c('blue', 'yellow'))(12), density.info='none',keysize=1)
dev.off()
heatmap.2(cor_matrix_mod, dendrogram='none', main ='all genes', Rowv=Null, Colv=Null, trace='none', col = colorRampPalette(c('blue', 'yellow'))(12), labRow=NA, labCol=NA , density.info='none',keysize=1)
write.table(cor_matrix_mod,"cor_matrix_mod.txt")

####
sample.markers <- FindAllMarkers(object = Fibroblast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
###
IB <- FindMarkers(Fibroblast, ident.1 = c("T3","T5","T10","T17","T18","T19","T21","T22","T24"), ident.2 = c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11"), verbose = FALSE, logfc.threshold=0, min.pct=0, min.cells.feature =1, min.cells.group =1)
write.table(IB,"IB_fibroblat_DE_all.txt",sep="\t")
IIA <- FindMarkers(Fibroblast, ident.1 = c("T6","T9"), ident.2 = c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11"), verbose = FALSE, logfc.threshold=0, min.pct=0, min.cells.feature =1, min.cells.group =1)
write.table(IIA,"IIA_fibroblat_DE_all.txt",sep="\t")
IIB <- FindMarkers(Fibroblast, ident.1 = c("T2","T4","T7","T11","T13","T14","T15","T16","T23"), ident.2 = c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11"), verbose = FALSE, logfc.threshold=0, min.pct=0, min.cells.feature =1, min.cells.group =1)
write.table(IIB,"IIB_fibroblat_DE_all.txt",sep="\t")
III <- FindMarkers(Fibroblast, ident.1 = c("T1","T8","T12"), ident.2 = c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11"), verbose = FALSE, logfc.threshold=0, min.pct=0, min.cells.feature =1, min.cells.group =1)
write.table(III,"III_fibroblat_DE_all.txt",sep="\t")
common_up_ELI <-intersect(intersect(intersect(intersect(up100gs, rownames(IB)), rownames(IIA)), rownames(IIB)), rownames(III))
common_dn_ELI <- intersect(intersect(intersect(intersect(dn100gs, rownames(IB)), rownames(IIA)), rownames(IIB)), rownames(III))
common_up <-intersect(intersect(intersect(rownames(IB), rownames(IIA)), rownames(IIB)), rownames(III))
common_dn <- intersect(intersect(intersect(rownames(IB), rownames(IIA)), rownames(IIB)), rownames(III))

merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}
common_up_ELI_merge <- merge.all(IB[common_up_ELI,], IIA[common_up_ELI,],IIB[common_up_ELI,],III[common_up_ELI,])
common_up_ELI_merge2<-common_up_ELI_merge[,c(2,7,12,17)]
colnames(common_up_ELI_merge2)<-c("IB","IIA","IIB","III")
common_up_ELI_merge2_sort_IB<-common_up_ELI_merge2[order(common_up_ELI_merge2$IB),]
write.table(common_up_ELI_merge2_sort_IB,"common_up_ELI_merge2_sort_IB.txt",sep="\t")

common_dn_ELI_merge <- merge.all(IB[common_dn_ELI,], IIA[common_dn_ELI,],IIB[common_dn_ELI,],III[common_dn_ELI,])
common_dn_ELI_merge2<-common_dn_ELI_merge[,c(2,7,12,17)]
colnames(common_dn_ELI_merge2)<-c("IB","IIA","IIB","III")
common_dn_ELI_merge2_sort_IB<-common_dn_ELI_merge2[order(common_dn_ELI_merge2$IB),]
write.table(common_dn_ELI_merge2_sort_IB,"common_dn_ELI_merge2_sort_IB.txt",sep="\t")


common_up_merge <- merge.all(IB[common_up,], IIA[common_up,],IIB[common_up,],III[common_up,])
common_up_merge2<-common_up_merge[,c(2,7,12,17)]
colnames(common_up_merge2)<-c("IB","IIA","IIB","III")
common_up_merge2_sort_IB<-common_up_merge2[order(common_up_merge2$IB),]
write.table(common_up_merge2_sort_IB,"common_up_merge2_sort_IB.txt",sep="\t")

common_dn_merge <- merge.all(IB[common_dn,], IIA[common_dn,],IIB[common_dn,],III[common_dn,])
common_dn_merge2<-common_dn_merge[,c(2,7,12,17)]
colnames(common_dn_merge2)<-c("IB","IIA","IIB","III")
common_dn_merge2_sort_IB<-common_dn_merge2[order(common_dn_merge2$IB),]
write.table(common_dn_merge2_sort_IB,"common_dn_merge2_sort_IB.txt",sep="\t")
###
common_dn_ELI_merge2_sort_IB <- as.matrix(read.table("common_dn_ELI_merge2_sort_IB.txt",header=T))
#heatmap.2(common_dn_ELI_merge2_sort_IB, dendrogram='none',main ='dn_100', Rowv=F, Colv=F, trace='none', col = colorRampPalette(c('blue', 'yellow'))(12), density.info='none',keysize=1)
common_dn_ELI_merge2_sort_IB <-as.data.frame(common_dn_ELI_merge2_sort_IB)
common_dn_ELI_merge2_sort_IB_add<-cbind(seq(1,dim(common_dn_ELI_merge2_sort_IB)[1]),common_dn_ELI_merge2_sort_IB)
colnames(common_dn_ELI_merge2_sort_IB_add)[1]<-c("geneorder")
common_dn_ELI_merge2_sort_IB_add<-as.data.frame(common_dn_ELI_merge2_sort_IB_add)
plot(common_dn_ELI_merge2_sort_IB_add$geneorder, common_dn_ELI_merge2_sort_IB_add$IB, type = "l", col = 1, xlab="genes", ylab="log2fc")
lines(common_dn_ELI_merge2_sort_IB_add$geneorder, common_dn_ELI_merge2_sort_IB_add$IIA, type = "l", col = 2)
lines(common_dn_ELI_merge2_sort_IB_add$geneorder, common_dn_ELI_merge2_sort_IB_add$IIB, type = "l", col = 3)
lines(common_dn_ELI_merge2_sort_IB_add$geneorder, common_dn_ELI_merge2_sort_IB_add$III, type = "l", col = 4)
legend(80,-0.4,legend=c("IB", "IIA","IIB","III"),col=c(1,2,3,4), lty=1,cex=0.8)

common_up_ELI_merge2_sort_IB <- as.matrix(read.table("common_up_ELI_merge2_sort_IB.txt",header=T))
#heatmap.2(common_up_ELI_merge2_sort_IB, dendrogram='none',main ='up_100', Rowv=F, Colv=F, trace='none', col = colorRampPalette(c('blue', 'yellow'))(12), density.info='none',keysize=1)
common_up_ELI_merge2_sort_IB <-as.data.frame(common_up_ELI_merge2_sort_IB)
common_up_ELI_merge2_sort_IB_add<-cbind(seq(1,dim(common_up_ELI_merge2_sort_IB)[1]),common_up_ELI_merge2_sort_IB)
colnames(common_up_ELI_merge2_sort_IB_add)[1]<-c("geneorder")
common_up_ELI_merge2_sort_IB_add<-as.data.frame(common_up_ELI_merge2_sort_IB_add)
plot(common_up_ELI_merge2_sort_IB_add$geneorder, common_up_ELI_merge2_sort_IB_add$IB, type = "l", col = 1, xlab="genes", ylab="log2fc")
lines(common_up_ELI_merge2_sort_IB_add$geneorder, common_up_ELI_merge2_sort_IB_add$IIA, type = "l", col = 2)
lines(common_up_ELI_merge2_sort_IB_add$geneorder, common_up_ELI_merge2_sort_IB_add$IIB, type = "l", col = 3)
lines(common_up_ELI_merge2_sort_IB_add$geneorder, common_up_ELI_merge2_sort_IB_add$III, type = "l", col = 4)
legend(80,-0.4,legend=c("IB", "IIA","IIB","III"),col=c(1,2,3,4), lty=1,cex=0.8)


common_dn_merge2_sort_IB <- as.matrix(read.table("common_dn_merge2_sort_IB.txt",header=T))
#heatmap.2(common_dn_merge2_sort_IB, dendrogram='none',main ='dn_100', Rowv=F, Colv=F, trace='none', col = colorRampPalette(c('blue', 'yellow'))(12), density.info='none',keysize=1)
common_dn_merge2_sort_IB <-as.data.frame(common_dn_merge2_sort_IB)
common_dn_merge2_sort_IB_add<-cbind(seq(1,dim(common_dn_merge2_sort_IB)[1]),common_dn_merge2_sort_IB)
colnames(common_dn_merge2_sort_IB_add)[1]<-c("geneorder")
common_dn_merge2_sort_IB_add<-as.data.frame(common_dn_merge2_sort_IB_add)
plot(common_dn_merge2_sort_IB_add$geneorder, common_dn_merge2_sort_IB_add$IB, type = "l", col = 1, xlab="genes", ylab="log2fc")
lines(common_dn_merge2_sort_IB_add$geneorder, common_dn_merge2_sort_IB_add$IIA, type = "l", col = 2)
lines(common_dn_merge2_sort_IB_add$geneorder, common_dn_merge2_sort_IB_add$IIB, type = "l", col = 3)
lines(common_dn_merge2_sort_IB_add$geneorder, common_dn_merge2_sort_IB_add$III, type = "l", col = 4)
legend(2000,-2,legend=c("IB", "IIA","IIB","III"),col=c(1,2,3,4), lty=1,cex=0.8)
###
common_dn_merge2<-common_dn_merge2_sort_IB
all_IB<-common_dn_merge2$IB
names(all_IB)<-rownames(common_dn_merge2)
all_IB_sort = sort(all_IB, decreasing = TRUE)
up100_IB<-intersect(up100gs, names(all_IB))
list_up100_IB<-list(up100_IB)
names(list_up100_IB)<-c("up100_IB")
fgsea_up <- fgsea(list_up100_IB,all_IB_sort,nperm=10000)
plotEnrichment(list_up100_IB[["up100_IB"]],all_IB_sort) + labs(title="up 100_IB")

all_IIA<-common_dn_merge2$IIA
names(all_IIA)<-rownames(common_dn_merge2)
all_IIA_sort = sort(all_IIA, decreasing = TRUE)
up100_IIA<-intersect(up100gs, names(all_IIA))
list_up100_IIA<-list(up100_IIA)
names(list_up100_IIA)<-c("up100_IIA")
fgsea_up <- fgsea(list_up100_IIA,all_IIA_sort,nperm=10000)
plotEnrichment(list_up100_IIA[["up100_IIA"]],all_IIA_sort) + labs(title="up 100_IIA")

all_IIB<-common_dn_merge2$IIB
names(all_IIB)<-rownames(common_dn_merge2)
all_IIB_sort = sort(all_IIB, decreasing = TRUE)
up100_IIB<-intersect(up100gs, names(all_IIB))
list_up100_IIB<-list(up100_IIB)
names(list_up100_IIB)<-c("up100_IIB")
fgsea_up <- fgsea(list_up100_IIB,all_IIB_sort,nperm=10000)
plotEnrichment(list_up100_IIB[["up100_IIB"]],all_IIB_sort) + labs(title="up 100_IIB")

all_III<-common_dn_merge2$III
names(all_III)<-rownames(common_dn_merge2)
all_III_sort = sort(all_III, decreasing = TRUE)
up100_III<-intersect(up100gs, names(all_III))
list_up100_III<-list(up100_III)
names(list_up100_III)<-c("up100_III")
fgsea_up <- fgsea(list_up100_III,all_III_sort,nperm=10000)
plotEnrichment(list_up100_III[["up100_III"]],all_III_sort) + labs(title="up 100_III")

###
common_dn_merge2<-common_dn_merge2_sort_IB
all_IB<-common_dn_merge2$IB
names(all_IB)<-rownames(common_dn_merge2)
all_IB_sort = sort(all_IB, decreasing = TRUE)
dn100_IB<-intersect(dn100gs, names(all_IB))
list_dn100_IB<-list(dn100_IB)
names(list_dn100_IB)<-c("dn100_IB")
fgsea_dn <- fgsea(list_dn100_IB,all_IB_sort,nperm=10000)
plotEnrichment(list_dn100_IB[["dn100_IB"]],all_IB_sort) + labs(title="dn 100_IB")

all_IIA<-common_dn_merge2$IIA
names(all_IIA)<-rownames(common_dn_merge2)
all_IIA_sort = sort(all_IIA, decreasing = TRUE)
dn100_IIA<-intersect(dn100gs, names(all_IIA))
list_dn100_IIA<-list(dn100_IIA)
names(list_dn100_IIA)<-c("dn100_IIA")
fgsea_dn <- fgsea(list_dn100_IIA,all_IIA_sort,nperm=10000)
plotEnrichment(list_dn100_IIA[["dn100_IIA"]],all_IIA_sort) + labs(title="dn 100_IIA")

all_IIB<-common_dn_merge2$IIB
names(all_IIB)<-rownames(common_dn_merge2)
all_IIB_sort = sort(all_IIB, decreasing = TRUE)
dn100_IIB<-intersect(dn100gs, names(all_IIB))
list_dn100_IIB<-list(dn100_IIB)
names(list_dn100_IIB)<-c("dn100_IIB")
fgsea_dn <- fgsea(list_dn100_IIB,all_IIB_sort,nperm=10000)
plotEnrichment(list_dn100_IIB[["dn100_IIB"]],all_IIB_sort) + labs(title="dn 100_IIB")

all_III<-common_dn_merge2$III
names(all_III)<-rownames(common_dn_merge2)
all_III_sort = sort(all_III, decreasing = TRUE)
dn100_III<-intersect(dn100gs, names(all_III))
list_dn100_III<-list(dn100_III)
names(list_dn100_III)<-c("dn100_III")
fgsea_dn <- fgsea(list_dn100_III,all_III_sort,nperm=10000)
plotEnrichment(list_dn100_III[["dn100_III"]],all_III_sort) + labs(title="dn 100_III")
#############
all_fc<-abs(ALL$avg_log2FC)
names(all_fc)<-rownames(ALL)
all_fc_sort = sort(all_fc, decreasing = TRUE)
dn100_fc<-intersect(dn100gs, names(all_fc))
list_dn100_fc<-list(dn100_fc)
names(list_dn100_fc)<-c("dn100_fc")
fgsea_dn <- fgsea(list_dn100_fc,all_fc_sort,nperm=10000)
plotEnrichment(list_dn100_fc[["dn100_fc"]],all_fc_sort) + labs(title="dn 100_fc")

all_fc<-ALL$avg_log2FC 
names(all_fc)<-rownames(ALL)
all_fc_sort = sort(all_fc, decreasing = TRUE)
up100_fc<-intersect(up100gs, names(all_fc))
list_up100_fc<-list(up100_fc)
names(list_up100_fc)<-c("up100_fc")
fgsea_up <- fgsea(list_up100_fc,all_fc_sort,nperm=10000)
plotEnrichment(list_up100_fc[["up100_fc"]],all_fc_sort) + labs(title="up 100_fc")
###cor-->from line 154
cor_matrix_mod<-cor(t(matrix_mod))
rest_of_genes = setdiff(setdiff(allgenes, common_up), common_dn)
write.table(cor_matrix_mod,"cor_matrix_mod.txt")
t.test(c(cor_matrix_mod[common_dn,common_dn]), c(cor_matrix_mod[rest_of_genes,rest_of_genes]))
t.test(c(cor_matrix_mod[common_up,common_up]), c(cor_matrix_mod[rest_of_genes,rest_of_genes]))

#################
Fibroblast <- subset(bothdata, idents = c("Fibroblast_TUMOR","Fibroblast_CTRL"))
Idents(object = Fibroblast) <- Fibroblast@meta.data$orig.ident
IIA_vs_IB <- FindMarkers(Fibroblast, ident.1 = c("T6","T9"), ident.2 = c("T3","T5","T10","T17","T18","T19","T21","T22","T24"), verbose = FALSE, logfc.threshold=0, min.pct=0, min.cells.feature =1, min.cells.group =1)
write.table(IIA_vs_IB,"IIA_vs_IB_fibroblat_DE_all.txt",sep="\t")
IIB_vs_IIA <- FindMarkers(Fibroblast, ident.1 = c("T2","T4","T7","T11","T13","T14","T15","T16","T23"), ident.2 = c("T6","T9"), verbose = FALSE, logfc.threshold=0, min.pct=0, min.cells.feature =1, min.cells.group =1)
write.table(IIB_vs_IIA,"IIB_vs_IIA_fibroblat_DE_all.txt",sep="\t")
III_vs_IIB <- FindMarkers(Fibroblast, ident.1 = c("T1","T8","T12"), ident.2 = c("T2","T4","T7","T11","T13","T14","T15","T16","T23"), verbose = FALSE, logfc.threshold=0, min.pct=0, min.cells.feature =1, min.cells.group =1)
write.table(III_vs_IIB,"III_vs_IIB_fibroblat_DE_all.txt",sep="\t")
III_vs_II <- FindMarkers(Fibroblast, ident.1 = c("T1","T8","T12"), ident.2 = c("T2","T4","T7","T11","T13","T14","T15","T16","T23","T6","T9"), verbose = FALSE, logfc.threshold=0, min.pct=0, min.cells.feature =1, min.cells.group =1)
write.table(III_vs_II,"III_vs_II_fibroblat_DE_all.txt",sep="\t")
#####
IIA_vs_IB<-read.table("IIA_vs_IB_fibroblat_DE_all.txt",header=TRUE)
list_up100gs<-list(up100gs)
names(list_up100gs)<-c("up100")
a<-as.vector(IIA_vs_IB$avg_log2FC)
names(a)<-rownames(IIA_vs_IB)
mcov_sort = sort(a, decreasing = TRUE)
mcov_sort = mcov_sort[!duplicated(names(mcov_sort))]
mcov_sort<-mcov_sort[!is.na(mcov_sort) & !is.infinite(mcov_sort)]
fgsea_up <- fgsea(list_up100gs,mcov_sort,nperm=10000)
plotEnrichment(list_up100gs[["up100"]],mcov_sort) + labs(title="up 100")

list_dn100gs<-list(dn100gs)
names(list_dn100gs)<-c("dn100")
fgsea_dn <- fgsea(list_dn100gs,mcov_sort,nperm=10000)
plotEnrichment(list_dn100gs[["dn100"]],mcov_sort) + labs(title="dn 100")
dev.off()


IIB_vs_IIA<-read.table("IIB_vs_IIA_fibroblat_DE_all.txt",header=TRUE)
list_up100gs<-list(up100gs)
names(list_up100gs)<-c("up100")
a<-as.vector(IIB_vs_IIA$avg_log2FC)
names(a)<-rownames(IIB_vs_IIA)
mcov_sort = sort(a, decreasing = TRUE)
mcov_sort = mcov_sort[!duplicated(names(mcov_sort))]
mcov_sort<-mcov_sort[!is.na(mcov_sort) & !is.infinite(mcov_sort)]
fgsea_up <- fgsea(list_up100gs,mcov_sort,nperm=10000)
plotEnrichment(list_up100gs[["up100"]],mcov_sort) + labs(title="up 100")

list_dn100gs<-list(dn100gs)
names(list_dn100gs)<-c("dn100")
fgsea_dn <- fgsea(list_dn100gs,mcov_sort,nperm=10000)
plotEnrichment(list_dn100gs[["dn100"]],mcov_sort) + labs(title="dn 100")
dev.off()

III_vs_II<-read.table("III_vs_II_fibroblat_DE_all.txt",header=TRUE)
list_up100gs<-list(up100gs)
names(list_up100gs)<-c("up100")
a<-as.vector(III_vs_II$avg_log2FC)
names(a)<-rownames(III_vs_II)
mcov_sort = sort(a, decreasing = TRUE)
mcov_sort = mcov_sort[!duplicated(names(mcov_sort))]
mcov_sort<-mcov_sort[!is.na(mcov_sort) & !is.infinite(mcov_sort)]
fgsea_up <- fgsea(list_up100gs,mcov_sort,nperm=10000)
plotEnrichment(list_up100gs[["up100"]],mcov_sort) + labs(title="up 100")

list_dn100gs<-list(dn100gs)
names(list_dn100gs)<-c("dn100")
fgsea_dn <- fgsea(list_dn100gs,mcov_sort,nperm=10000)
plotEnrichment(list_dn100gs[["dn100"]],mcov_sort) + labs(title="dn 100")
dev.off()

III_vs_IIB<-read.table("III_vs_IIB_fibroblat_DE_all.txt",header=TRUE)
list_up100gs<-list(up100gs)
names(list_up100gs)<-c("up100")
a<-as.vector(III_vs_IIB$avg_log2FC)
names(a)<-rownames(III_vs_IIB)
mcov_sort = sort(a, decreasing = TRUE)
mcov_sort = mcov_sort[!duplicated(names(mcov_sort))]
mcov_sort<-mcov_sort[!is.na(mcov_sort) & !is.infinite(mcov_sort)]
fgsea_up <- fgsea(list_up100gs,mcov_sort,nperm=10000)
plotEnrichment(list_up100gs[["up100"]],mcov_sort) + labs(title="up 100")

list_dn100gs<-list(dn100gs)
names(list_dn100gs)<-c("dn100")
fgsea_dn <- fgsea(list_dn100gs,mcov_sort,nperm=10000)
plotEnrichment(list_dn100gs[["dn100"]],mcov_sort) + labs(title="dn 100")
dev.off()
