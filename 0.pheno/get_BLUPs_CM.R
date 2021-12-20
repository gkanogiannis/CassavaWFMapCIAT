#!/home/agkanogiannis/bin/Rscript

#
# get_BLUPs_CM.R
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
BiocManager::install(update=F)
if (!require(gdsfmt)) BiocManager::install(c("gdsfmt")) 
library(gdsfmt)
if (!require(SNPRelate)) BiocManager::install(c("SNPRelate"))
library(SNPRelate)
library(dplyr)
library(car)
library(lme4)
library(lmerTest)

pheno.all <- data.frame(clone=character(),rep=character(),plant=character(),leaf=character(),stringsAsFactors=F)
experiments.names <- list("CM_cnt2013_Exp1_Ane","CM_cnt2016_Exp2_Ane","CM_cnt2017_Exp3_Ane","CM_cnt2017_Exp4_Ane",
                          "CM_cnt2018_Exp5_Ane","CM_are2017_Exp3_Ane","CM_are2017_Exp4_Ane","CM_are2018_Exp5_Ane")
for(experiment in experiments.names){
  pheno.experiment <- read.table(paste0(experiment,".tsv"),sep="\t",header=T,
                                 colClasses=c("NULL","NULL","NULL","NULL","NULL",
                                   "character","character","character","character","numeric"),stringsAsFactors=F)
  names(pheno.experiment)[names(pheno.experiment) == "nymphs_area"] <- experiment
  #pheno.experiment <- aggregate(as.formula(paste("nymphs_area ~ clone+rep+plant+leaf")),pheno.experiment, mean)
  pheno.all <- merge(pheno.all,pheno.experiment,all=T,by=c("clone","rep","plant","leaf"))
}
pheno.all <- pheno.all %>% distinct()
write.table(pheno.all,file="pheno_CM.all.tsv",sep="\t",quote=F,row.names=F)

BLUP.all <- data.frame(clone=character(),stringsAsFactors=F)
pdf(file = "BLUP_CM.plots.pdf",paper="usr",width=11, height=8.5)
for(experiment in experiments.names){
  #experiment <- "CM_cnt2013_Exp1_Ane"
  layout(matrix(c(1,2,3,4,8,5,6,7), 2, 4, byrow=T))
  #layout.show(5)
  pheno.experiment <- pheno.all[!is.na(pheno.all[experiment]),c("clone","rep","plant","leaf",experiment)]
  pheno.experiment <- pheno.experiment %>% distinct()
  pheno.experiment$clone <- as.factor(pheno.experiment$clone)
  pheno.experiment$rep <- as.factor(pheno.experiment$rep)
  pheno.experiment$plant <- as.factor(pheno.experiment$plant)
  pheno.experiment$leaf <- as.factor(pheno.experiment$leaf)
  
  #par(cex.axis=0.5)
  #par(cex.lab=0.75)
  hist1<-hist(pheno.experiment[[experiment]], col="gold",xlab=experiment, main=paste0(experiment," histogram"), breaks=100)
  bp1<-Boxplot(pheno.experiment[[experiment]]~pheno.experiment$rep,id=list(n=Inf,cex=0.1), xlab=NULL, ylab=NULL, main=paste0(experiment,"\nby rep"), col="pink")
  bp2<-Boxplot(pheno.experiment[[experiment]]~pheno.experiment$plant,id=list(n=Inf,cex=0.1), xlab=NULL, ylab=NULL, main=paste0(experiment,"\nby plant"), col="pink")
  bp3<-Boxplot(pheno.experiment[[experiment]]~pheno.experiment$leaf,id=list(n=Inf,cex=0.1), xlab=NULL, ylab=NULL, main=paste0(experiment,"\nby leaf"), col="pink")
  
  #pheno.experiment <- pheno.experiment[-as.numeric(bp1),]
  #pheno.experiment <- pheno.experiment[-as.numeric(bp2),]
  #pheno.experiment <- pheno.experiment[-as.numeric(bp3),]
  
  bp1<-Boxplot(pheno.experiment[[experiment]]~pheno.experiment$rep,id=list(n=Inf,cex=0.1), xlab=NULL, ylab=NULL, main=NULL, col="pink")
  bp2<-Boxplot(pheno.experiment[[experiment]]~pheno.experiment$plant,id=list(n=Inf,cex=0.1), xlab=NULL, ylab=NULL, main=NULL, col="pink")
  bp3<-Boxplot(pheno.experiment[[experiment]]~pheno.experiment$leaf,id=list(n=Inf,cex=0.1), xlab=NULL, ylab=NULL, main=NULL, col="pink")
  
  model <- lmer(as.formula(paste0(experiment,
    " ~ 1 + 
    (1|clone)+(1|rep)+(1|plant)+(1|leaf)+
    (1|clone:rep)+(1|clone:plant)+(1|clone:leaf)+
    (1|rep:plant)+(1|rep:leaf)+(1|plant:leaf)")),
    data=pheno.experiment)
  
  #print(summary(model))
  BLUP <- ranef(model,drop=F)
  BLUP.clone <- data.frame(clone=rownames(BLUP$clone),experiment=BLUP$clone,stringsAsFactors=F)  
  colnames(BLUP.clone) <- c("clone",experiment)
  BLUP.clone[[experiment]] <- BLUP.clone[[experiment]] + fixef(model)[[1]]
  
  write.table(BLUP.clone,file=paste0("BLUPs.",experiment,".csv"),sep="\t",quote=F,row.names=F,col.names=T)
  BLUP.all <- merge(BLUP.all,BLUP.clone,all=T,by=c("clone"))
  
  ## Creating plots with the BLUPs
  # Create a histogram with the BLUP for each clone
  hist(BLUP.clone[,2], col="brown",xlab=paste0(experiment," BLUP"),main=paste0(experiment," BLUP histogram"),breaks=100)
}
dev.off()
write.table(BLUP.all,file=paste0("BLUPs_CM.csv"),sep="\t",quote=F,row.names=F,col.names=T)
