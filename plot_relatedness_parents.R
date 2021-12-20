#!/home/agkanogiannis/bin/Rscript

#
# plot_relatedness_parents.R
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

input <- "CM8996.filtered.relatedness2"
parent1 <- "COL2246_MALE"
parent2 <- "ECU72_FEMALE"

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  input <- args[1]
  parent1 <- args[2]
  parent2 <- args[3]
}

input_type = summary( file(input) )$class
print(input_type)

if(input_type=='gzfile'){
  data <- read.table(skip=0,file=gzfile(input),header=T,stringsAsFactors=F,fill=T)
}else{
  data <- read.table(skip=0,file=input,header=T,stringsAsFactors=F,fill=T)
}

data_parent1<-data[data$INDV1==parent1,c(1,2,7)]
data_parent1<-data_parent1[order(data_parent1$INDV2),]
data_parent2<-data[data$INDV1==parent2,c(1,2,7)]
data_parent2<-data_parent2[order(data_parent2$INDV2),]

relatedness<-cbind(data_parent1$RELATEDNESS_PHI,data_parent2$RELATEDNESS_PHI)

colnames(relatedness)<-c(parent1,parent2)
rownames(relatedness)<-data_parent1[,2]

library(RColorBrewer)
cols <- brewer.pal(8, "Dark2")
pdf(file=paste0(input,".pdf"),paper="us",width=11.0, height=8.5)
par(mar=c(3,3,3,3)+2,pty="s")
plot(relatedness[,parent2] ~ relatedness[,parent1], main=basename(input),
     ylab = parent2, xlab=parent1,
     xlim=c(0,0.5),ylim=c(0,0.5), col=NULL)
text(x=relatedness[,parent1], y=relatedness[,parent2], labels=rownames(relatedness),
       col=cols[1], pos=3, cex=0.25, offset=0.2)
points(x=relatedness[,parent1], y=relatedness[,parent2],
       col=cols[4], pch=".")
#grid(lty=1)
abline(h=0.354, lty=2)
abline(v=0.354, lty=2)
abline(h=0.177, lty=4)
abline(v=0.177, lty=4)
abline(h=0.0884, lty=3)
abline(v=0.0884, lty=3)
abline(h=0.0442, lty=3)
abline(v=0.0442, lty=3)
dev.off()

