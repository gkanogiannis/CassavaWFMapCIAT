#!/home/agkanogiannis/bin/Rscript

#
# relatedness.R
# Script part of CassavaWFMapCIAT project
#
# Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

if (!requireNamespace("BiocManager")) install.packages("BiocManager")
#BiocManager::install()
if (!require(gdsfmt)) BiocManager::install(c("gdsfmt"))
suppressPackageStartupMessages(library(gdsfmt))
if (!require(SNPRelate)) BiocManager::install(c("SNPRelate"))
suppressPackageStartupMessages(library(SNPRelate))
if (!require(Ternary)) install.packages("Ternary")
suppressPackageStartupMessages(library(Ternary))
if (!require(dynamicTreeCut)) install.packages("dynamicTreeCut")
suppressPackageStartupMessages(library(dynamicTreeCut))
if (!require(ape)) install.packages("ape")
suppressPackageStartupMessages(library(ape))
if (!require(dendextend)) install.packages("dendextend")
suppressPackageStartupMessages(library(dendextend))
if (!require(ggplot2)) install.packages("ggplot2")
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(lattice))

c25 <- c("dodgerblue2","#E31A1C","green4","#6A3D9A","#FF7F00","black",
         "gold1","skyblue2","#FB9A99","palegreen2","#CAB2D6","#FDBF6F",
         "gray70","khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
         "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown"
)
palette(c25)

args = commandArgs(trailingOnly=TRUE)
# 1 vcf, 2 pihat, 3 ld, 4 mle, 5 pop

pihat.threshold <- as.numeric(args[2])
if(is.na(pihat.threshold) | pihat.threshold>=1){
  pihat.threshold <- 0.9
}
kinship.threshold <- pihat.threshold/2

ld.threshold <- as.numeric(args[3])
if(is.na(ld.threshold) | ld.threshold>=1){
  ld.threshold <- 1.0
}

ibd.mle.calc=F
if(is.na(as.logical(args[4]))){
  ibd.mle.calc=F
}else{
  ibd.mle.calc=as.logical(args[4])
}

if(!is.na(args[1])){
  vcf <- args[1]
  #vcf <- "Cassava_RAD_CIAT_GWAS.snps.vmiss_0.10.smiss_0.25.MAF_0.03.hardy_1e-10.het_clean.ld.vcf"
  gds.file <- paste0(vcf,".gds")
  report.file <- paste0(vcf,".relatedness.pdf")
  snpgdsVCF2GDS(vcf, gds.file, method="biallelic.only", ignore.chr.prefix="Chromosome")
  snpgdsSummary(gds.file)
  genofile <- snpgdsOpen(gds.file)
  sample.id <- read.gdsn(index.gdsn(genofile, path="sample.id"))
  snp.id <- read.gdsn(index.gdsn(genofile, path="snp.id"))
}else{
  vcf <- NA
}
if(!is.na(args[5])){
  #pop.groups <- read.table("samples.RAD_all.sample_group.txt",header=F,stringsAsFactors=F,col.names=c("Genotype","Group"))
  #pop_code <- pop.groups[match(sample.id, pop.groups$Genotype),"Group"]
  pop.groups <- read.table(args[5],header=F,stringsAsFactors=F,col.names=c("Genotype","Group"))
  pop_code <- pop.groups[match(sample.id, pop.groups$Genotype),"Group"]
}else{
  pop_code <- NA
}

print(ibd.mle.calc)
print(pop_code)

#cairo_pdf(file=report.file,width=11, height=8.5)
pdf(file=report.file,width=11, height=8.5,paper="usr")
#par(pty="s")


# LD-based SNP pruning
if(ld.threshold<1){
  set.seed(1000)
  snp.ld <- snpgdsLDpruning(genofile,autosome.only=F,ld.threshold=ld.threshold)
  snp.ld.id <- unlist(unname(snp.ld))
} else{
  snp.ld.id <- snp.id
}


# Relatedness Analysis

# Fun ction that creates clustering from dist and cuts the tree at the cutheight
# Branches bellow the cutheight are considered clones
FunctionGetClustersClones <- function (dist,cutHeight,fileprefix,title,pop_code) {
  par_temp <- par()
  if(!is.na(pop_code)) {
    layout(matrix(c(1,2,3), 3, 1, byrow=TRUE))
  }else {
    layout(matrix(c(1,2), 2, 1, byrow=TRUE))
  }
  # Get tree form distance matrix
  tree <- hclust(dist,method="complete")
  # Write tree as newick
  write.tree(as.phylo(tree),paste0(fileprefix,".tree"),digits=4)
  # Cut tree dynamic at the pihat threshold to get clusters
  # Cluster 0 contains singletons (unrelated)
  cutree <- cutreeDynamic(tree,distM=as.matrix(dist),method="hybrid",
                          cutHeight=cutHeight,minClusterSize=1)
  clusters<-data.frame(label=tree$labels,cluster=cutree)
  rownames(clusters)<-tree$labels
  # Plot clusters cut at the cutheigt (clones)
  # (uniques blue, clone clusters other colors)
  hc <- snpgdsHCluster(as.matrix(dist),need.mat=T,hang=0)
  rv <- snpgdsCutTree(hc,label.H=T,samp.group=as.factor(clusters$cluster),pch.list=c(20))
  print(table(as.factor(clusters$cluster)))
  rv$dendrogram %>% 
    set("labels_cex", 0.1) %>%
    set("branches_lty", 3) %>%
    set("branches_lwd", 0.2) %>%
    set("leaves_cex",0.1) %>%
    plot(main=paste0(title," Clustering (Clones)"),dLeaf=0.2)
  abline(h=cutHeight, col="red",lty=3)
  text(0.5,cutHeight*1.2, paste0("cutheight=",cutHeight),cex=0.5)
  # Write clusters as tsv
  write.table(clusters[order(clusters$cluster),],sep="\t",
              file=paste0(fileprefix,".clusters.clones.tsv"),quote=F,row.names=F)
  # Write clones as tsv
  output <- ""
  for(i in 1:max(clusters$cluster)){
    sub <- subset(clusters,cluster==i)
    sub.len <- nrow(sub)
    if(sub.len>1){
      for(j in 1:sub.len) output <- paste0(output,sub[j,"label"],"\t")
      output <- paste0(output,"\n")
    }
  }
  write(output,paste0(fileprefix,".clones.tsv"))
  # Plot clusters cut automatically
  cutree <- cutreeDynamic(tree,distM=as.matrix(dist),method="hybrid")
  clusters<-data.frame(label=tree$labels,cluster=cutree)
  rownames(clusters)<-tree$labels
  hc <- snpgdsHCluster(as.matrix(dist),need.mat=T,hang=0)
  rv <- snpgdsCutTree(hc,label.H=T,samp.group=as.factor(clusters$cluster))
  print(table(rv$samp.group))
  rv$dendrogram %>% 
    set("labels_cex", 0.1) %>%
    set("branches_lty", 3) %>%
    set("branches_lwd", 0.2) %>%
    set("leaves_cex",0.1) %>%
    plot(main=paste0(title," Clustering (Auto)"),dLeaf=0.2)
  legend("topright", legend=levels(as.factor(clusters$cluster)),
         cex=0.5,pt.cex=1, pch=20, col=1:nlevels(as.factor(clusters$cluster)), ncol=4)
  # Write clusters as tsv
  write.table(clusters[order(clusters$cluster),],sep="\t",
              file=paste0(fileprefix,".clusters.auto.tsv"),quote=F,row.names=F)
  # Plot clusters if we have pop information
  if(!is.na(pop_code)){
    hc <- snpgdsHCluster(as.matrix(dist),need.mat=T,hang=0)
    rv <- snpgdsCutTree(hc,label.H=T,samp.group=as.factor(pop_code))
    print(table(rv$samp.group))
    rv$dendrogram %>% 
      set("labels_cex", 0.1) %>%
      set("branches_lty", 3) %>%
      set("branches_lwd", 0.2) %>%
      set("leaves_cex",0.1) %>%
      plot(main=paste0(title," Clustering (Pop)"),dLeaf=0.2)
    legend("topright", legend=levels(as.factor(pop_code)),
           cex=0.5,pt.cex=1, pch=20, col=1:nlevels(as.factor(pop_code)), ncol=4)
  }
  par(par_temp)
}
FunctionToContour <- function (a, b, c) {0.5*a+b}

############################
# IBD Analysis using PLINK method of moments (MoM) and Maximum Likelihood Estimation (MLE)
ibd.mom <- snpgdsIBDMoM(genofile, sample.id=sample.id, snp.id=snp.ld.id,kinship=T,autosome.only=F)
#ibd.mom$kinship[is.na(ibd.mom$kinship)] <- 0
ibd.mom.coeff <- snpgdsIBDSelection(ibd.mom)
ibd.mom.coeff <- cbind(ibd.mom.coeff,data.frame(k2=1-ibd.mom.coeff$k0-ibd.mom.coeff$k1))
ibd.mom.coeff <- cbind(ibd.mom.coeff,data.frame(pihat=0.5*ibd.mom.coeff$k1+ibd.mom.coeff$k2))
#ibd.mom.coeff <- ibd.mom.coeff[!is.na(ibd.mom.coeff$k0) & !is.na(ibd.mom.coeff$k1) & !is.na(ibd.mom.coeff$k2),]
if(ibd.mle.calc){
  set.seed(100)
  snp.mle.id <- sample(snp.ld.id, 1500)  # random 1500 SNPs
  ibd.mle <- snpgdsIBDMLE(genofile, sample.id=sample.id, snp.id=snp.mle.id,kinship=T,autosome.only=F,num.thread=4)
  ibd.mle$kinship[is.na(ibd.mle$kinship)] <- 0
  ibd.mle.coeff <- snpgdsIBDSelection(ibd.mle)
  ibd.mle.coeff <- cbind(ibd.mle.coeff,data.frame(k2=1-ibd.mle.coeff$k0-ibd.mle.coeff$k1))
  ibd.mle.coeff <- cbind(ibd.mle.coeff,data.frame(pihat=0.5*ibd.mle.coeff$k1+ibd.mle.coeff$k2))
  ibd.mle.coeff <- ibd.mle.coeff[!is.na(ibd.mle.coeff$k0) & !is.na(ibd.mle.coeff$k1) & !is.na(ibd.mle.coeff$k2),]
}
# Plot Cotterman coeffs
par_temp <- par()
if(ibd.mle.calc) {
  layout(matrix(c(1,2), 1, 2, byrow=TRUE))
}else {
  layout(matrix(c(1), 1, 1, byrow=TRUE))
}
TernaryPlot(point="up",lab.cex=0.5, grid.minor.lines=0,
            grid.lty='solid', col=rgb(0.9, 0.9, 0.9), grid.col="white", 
            axis.col=rgb(0.6, 0.6, 0.6), ticks.col=rgb(0.6, 0.6, 0.6),
            padding=0.08,main="IBD (MoM) Cotterman",
            alab="k1 \u2192", blab="k2 \u2192", clab="\u2190 k0",
            atip="k1", btip="k2", ctip="k0",clockwise=T)
ColourTernary(TernaryPointValues(FunctionToContour, resolution=24L))
TernaryContour(FunctionToContour, resolution=24L)
AddToTernary(points,ibd.mom.coeff[,c("k1","k2","k0")],pch=".", cex=1.0)
if(ibd.mle.calc){
  TernaryPlot(point="up",lab.cex=0.5, grid.minor.lines=0,
              grid.lty='solid', col=rgb(0.9, 0.9, 0.9), grid.col="white", 
              axis.col=rgb(0.6, 0.6, 0.6), ticks.col=rgb(0.6, 0.6, 0.6),
              padding=0.08,main="IBD (MLE) Cotterman",
              alab="k1 \u2192", blab="k2 \u2192", clab="\u2190 k0",
              atip="k1", btip="k2", ctip="k0",clockwise=T)
  ColourTernary(TernaryPointValues(FunctionToContour, resolution=24L))
  TernaryContour(FunctionToContour, resolution=24L)
  AddToTernary(points,ibd.mle.coeff[,c("k1","k2","k0")],pch=".", cex=1.0)
}
# Get distance and set kinship to zero if negative
ibd.mom$kinship[ibd.mom$kinship < 0] <- 0
ibd.mom.coeff$kinship[ibd.mom.coeff$kinship<0] <- 0
ibd.mom.dist <- 1 - 2.0*ibd.mom$kinship
colnames(ibd.mom.dist)<-ibd.mom$sample.id
rownames(ibd.mom.dist)<-ibd.mom$sample.id
# Write dist
write(length(ibd.mom$sample.id), file=paste0(vcf,".ibd_mom.dist"),sep="\n")
write.table(ibd.mom.dist,file=paste0(vcf,".ibd_mom.dist"),quote=F,sep="\t",col.names=F,append=T)
if(ibd.mle.calc){
  ibd.mle$kinship[ibd.mle$kinship < 0] <- 0
  ibd.mle.coeff$kinship[ibd.mle.coeff$kinship<0] <- 0
  ibd.mle.dist <- 1 - 2.0*ibd.mle$kinship
  colnames(ibd.mle.dist)<-ibd.mle$sample.id
  rownames(ibd.mle.dist)<-ibd.mle$sample.id
  # Write dist
  write(length(ibd.mle$sample.id), file=paste0(vcf,".ibd_mle.dist"),sep="\n")
  write.table(ibd.mle.dist,file=paste0(vcf,".ibd_mle.dist"),quote=F,sep="\t",col.names=F,append=T)
}
# Get clusters and clones
FunctionGetClustersClones(as.dist(ibd.mom.dist),cutHeight=(1-pihat.threshold),
                          file=paste0(vcf,".ibd_mom"),title="IBD MoM",pop_code=pop_code)
if(ibd.mle.calc){
  FunctionGetClustersClones(as.dist(ibd.mle.dist),cutHeight=(1-pihat.threshold),
                            file=paste0(vcf,".ibd_mle"),title="IBD MLE",pop_code=pop_code)
}
# Write close related as tsv
samples.close_rel <- ibd.mom.coeff[ibd.mom.coeff$pihat>=pihat.threshold,c("ID1","ID2","kinship","pihat")]
write.table(samples.close_rel,file=paste0(vcf,".ibd_mom.close_related.tsv"),quote=F,sep="\t",row.names=F)
if(ibd.mle.calc){
  samples.close_rel <- ibd.mle.coeff[ibd.mle.coeff$pihat>=pihat.threshold,c("ID1","ID2","kinship","pihat")]
  write.table(samples.close_rel,file=paste0(vcf,".ibd_mle.close_related.tsv"),quote=F,sep="\t",row.names=F)
}

############################
# IBS Analysis
ibs <- snpgdsIBS(genofile,sample.id=sample.id,snp.id=snp.ld.id,autosome.only=F)
# Get distance
ibs.dist <- 1 - (ibs$ibs-min(ibs$ibs))/(max(ibs$ibs)-min(ibs$ibs))
colnames(ibs.dist)<-ibs$sample.id
rownames(ibs.dist)<-ibs$sample.id
# Write dist
write(length(ibs$sample.id), file=paste0(vcf,".ibs.dist"),sep="\n")
write.table(ibs.dist,file=paste0(vcf,".ibs.dist"),quote=F,sep="\t",col.names=F,append=T)
# Get clusters and clones
FunctionGetClustersClones(as.dist(ibs.dist),cutHeight=(1-pihat.threshold),
                          file=paste0(vcf,".ibs"),title="IBS",pop_code=pop_code)
# Write close related as tsv
samples.close_rel <- "ID1\tID2\tIBS\n"
for(i in ibs$sample.id)
  for(j in ibs$sample.id)
    if(i!=j & ibs.dist[i,j]<(1-pihat.threshold))
      samples.close_rel <- paste0(samples.close_rel,(paste0(i,"\t",j,"\t",1-ibs.dist[i,j],"\n")))
write(samples.close_rel,file=paste0(vcf,".ibs.close_related.tsv"))

############################
# PCA and Multidimensional Scaling Analysis (IBS Distance) plots
pca <- snpgdsPCA(genofile,autosome.only=F,snp.id=snp.ld.id,algorithm="exact",eigen.method="DSPEV")
# Perform multidimensional scaling analysis on the n×n matrix of genome-wide IBS pairwise distances:
ibs.loc <- cmdscale(ibs.dist, k = 2)
par_temp <- par()
layout(matrix(c(1,2), 1, 2, byrow=TRUE))
if(is.na(pop_code)){
  # In there is no population information
  plot.pca <- 
    ggplot(data.frame(EV1=pca$eigenvect[,1],EV2=pca$eigenvect[,2]),aes(EV1, EV2)) + coord_fixed() +
    geom_point(size=0.1) + 
    geom_text(aes(label=pca$sample.id),hjust=-0.1,vjust=0,size=0.1) +
    labs(title="PCA", x="eigenvector 1",y="eigenvector 2")
  plot.ibs <- 
    ggplot(data.frame(x=ibs.loc[, 1],y=ibs.loc[, 2]), aes(x, y)) + coord_fixed() +
    geom_point(size=0.1) +
    geom_text(aes(label=rownames(ibs.loc)),hjust=-0.1,vjust=0,size=0.1) +
    labs(title="IBS \n(Multidimensional Scaling Analysis)",x="",y="")
  print(plot.pca)
  print(plot.ibs)
}else{
  # If there is population information,
  # assume the order of sample IDs is as the same as population codes
  plot.pca_pop <- 
    ggplot(data.frame(EV1=pca$eigenvect[,1],EV2=pca$eigenvect[,2]),aes(EV1, EV2)) + coord_fixed() +
    geom_point(aes(color = factor(pop_code)),size=0.1) +
    geom_text(aes(label=pca$sample.id),hjust=-0.1,vjust=0,size=0.1) +
    labs(title="PCA", x="eigenvector 1",y="eigenvector 2",color = "Pop") + 
    theme(legend.justification = c("left", "top")) + scale_colour_manual(values=c25)
  plot.ibs_pop <- 
    ggplot(data.frame(x=ibs.loc[, 1],y=ibs.loc[, 2]), aes(x, y)) + coord_fixed() +
    geom_point(aes(color = factor(pop_code)),size=0.1) +
    geom_text(aes(label=rownames(ibs.loc)),hjust=-0.1,vjust=0,size=0.1) +
    labs(title="IBS \n(Multidimensional Scaling Analysis with pop)",x="",y="",color = "Pop") + 
    theme(legend.justification = c("left", "top")) + scale_colour_manual(values=c25)
  print(plot.pca_pop)
  print(plot.ibs_pop)
}
par(par_temp)


############################
# Estimating IBD Using KING method
#ibd.king <- snpgdsIBDKING(genofile,sample.id=sample.id,snp.id=snp.ld.id,autosome.only=F,type="KING-robust")
#ibd.king.coeff <- snpgdsIBDSelection(ibd.king)
## Get distance and set kinship to zero if negative
#ibd.king$kinship[ibd.king$kinship < 0] <- 0
#ibd.king.coeff$kinship[ibd.king.coeff$kinship<0] <- 0
#ibd.king.dist <- 1 - 2.0*ibd.king$kinship
#colnames(ibd.king.dist)<-ibd.king$sample.id
#rownames(ibd.king.dist)<-ibd.king$sample.id
## Plot Multidimensional Scaling Analysis (IBD kinship)
#par_temp <- par()
#if(!is.na(pop_code)) {
#  layout(matrix(c(1,2), 1, 2, byrow=TRUE))
#}else {
#  layout(matrix(c(1), 1, 1, byrow=TRUE))
#}
#loc <- cmdscale(ibd.king.dist, k = 2)
#x <- loc[, 1]; y <- loc[, 2]
#plot(x, y, xlab = "", ylab = "", pch=20, cex=1.0,
#     main = "IBD KING (Multidimensional Scaling Analysis)")
## MSA with pop
#if(!is.na(pop_code)){
#  # To perform multidimensional scaling analysis on the n×n matrix of genome-wide IBS pairwise distances:
#  loc <- cmdscale(ibd.king.dist, k = 2)
#  x <- loc[, 1]; y <- loc[, 2]
#  plot(x, y, col=as.integer(as.factor(pop_code)), xlab = "", ylab = "", pch=20, cex=1.0,
#       main = "IBD KING (Multidimensional Scaling Analysis with pop)")
#  legend("bottomleft", legend=levels(as.factor(pop_code)),cex=0.5,pt.cex=1, pch=20, col=1:nlevels(as.factor(pop_code)))
#}
#par(par_temp)
## Get clusters and clones
#FunctionGetClustersClones(as.dist(ibd.king.dist),cutHeight=(1-pihat.threshold),
#                          file=paste0(vcf,".ibd_king"),title="IBD KING",pop_code=pop_code)
## Write close related as tsv
#samples.close_rel <- ibd.king.coeff[ibd.king.coeff$kinship>=kinship.threshold,c("ID1","ID2","kinship")]
#write.table(samples.close_rel,file=paste0(vcf,".ibd_king.close_related.tsv"),quote=F,sep="\t",row.names=F)


dev.off()

snpgdsClose(genofile)
if (file.exists(gds.file)) file.remove(gds.file)
