#!/home/agkanogiannis/bin/Rscript

#
# plot_dist_parents.R
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

dist_input <- "CM8996.filtered.dist"
parent1 <- "COL2246_MALE"
parent2 <- "ECU72_FEMALE"

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  dist_input <- args[1]
  parent1 <- args[2]
  parent2 <- args[3]
}

dist_input_type = summary( file(dist_input) )$class
print(dist_input_type)

if(dist_input_type=='gzfile'){
  data <- read.table(skip=1,file=gzfile(dist_input),header=F,stringsAsFactors=F,fill=T)
}else{
  data <- read.table(skip=1,file=dist_input,header=F,stringsAsFactors=F,fill=T)
}

rownames(data) <- data[,1]
data[,1] <- NULL
colnames(data) <- rownames(data)

dist<-cbind(data[parent1],data[parent2])

library(RColorBrewer)
cols <- brewer.pal(8, "Dark2")
pdf(file=paste0(dist_input,".pdf"),paper="us",width=11.0, height=8.5)
par(mar=c(3,3,3,3)+2,pty="s")
plot(dist[,parent2] ~ dist[,parent1], main=basename(dist_input),
     ylab = parent2, xlab=parent1,
     xlim=c(0,1.0),ylim=c(0,1.0), col=NULL)
text(x=dist[,parent1], y=dist[,parent2], labels=rownames(dist),
       col=cols[1], pos=3, cex=0.25, offset=0.2)
points(x=dist[,parent1], y=dist[,parent2],
       col=cols[4], pch=".")
#grid(lty=1)
abline(h=0.65, lty=2)
abline(v=0.66, lty=2)
abline(h=0.72, lty=4)
abline(v=0.72, lty=4)
abline(h=0.79, lty=3)
abline(v=0.79, lty=3)
dev.off()

