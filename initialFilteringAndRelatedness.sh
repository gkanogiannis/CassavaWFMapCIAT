#!/bin/bash

#
# initialFilteringAndRelatedness.sh
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

usage() { echo "Usage: $0 [-i <input_file>] [-v <vmiss>] [-s <smiss>] [-m <maf MIN>] [-M <maf MAX>] [-h <hwe>] [-f <vcf INFO stats>] [-p <pihat>] [-l <ld>] [-g <pop_groups] [-e <use MLE>]" 1>&2; exit 1; }

vmiss=1 ; #0.25 
smiss=1 ; #0.5
maf=0 ; #0.01
maf_max=0.5 ;
hwe=0 ; #1e-10
hwe_tmp=0 ;
info_stats=0 ;
pihat=1 ; #0.9
ld=1 ; #0.2
pop_groups="" ;
MLE="false" ;

while getopts ":i:v:s:m:M:h:fp:l:g:e" opt; do
	case "${opt}" in
		i)
			input_file=${OPTARG}
			;;
		v)
			vmiss=${OPTARG}
			(( `echo "${vmiss}>=0 && ${vmiss}<=1"|bc -l` )) || usage
			;;
		s)
			smiss=${OPTARG}
			(( `echo "${smiss}>=0 && ${smiss}<=1"|bc -l` )) || usage
			;;
		m)
			maf=${OPTARG}
			(( `echo "${maf}>=0 && ${maf}<=0.5"|bc -l` )) || usage
			;;
		M)
			maf_max=${OPTARG}
			(( `echo "${maf_max}>=0 && ${maf_max}<=0.5"|bc -l` )) || usage
			;;
		h)
			hwe=${OPTARG}
			hwe_tmp=`sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"${hwe}"`			
			(( `echo "${hwe_tmp}>=0 && ${hwe_tmp}<=1"|bc -l` )) || usage
			;;
		f)
			info_stats=1
			;;
		p)
			pihat=${OPTARG}
			(( `echo "${pihat}>=0 && ${pihat}<=1"|bc -l` )) || usage
			;;
		l)
			ld=${OPTARG}
			(( `echo "${ld}>=0 && ${ld}<=1"|bc -l` )) || usage
			;;
		g)
			pop_groups=${OPTARG}
			;;
		e)
			MLE="true"
			;;
		*)
			usage
			;;
	esac
done
shift $((OPTIND-1))


echo -e "input_file=${input_file}\nvmiss=${vmiss}\nsmiss=${smiss}\nmaf=${maf}\nhwe=${hwe}\ninfo_stats=${info_stats}\npihat=${pihat}\nld=${ld}\npop_groups=${pop_groups}\nMLE=${MLE}"

unset R_HOME ;

#basename=Cassava_RAD_CIAT_GWAS.snps
basename=`echo -e "${input_file}" | sed 's/.vcf.*//'` ;


relatedness_R=/home/agkanogiannis/CIAT/scripts/relatedness.R ;
vcftools=/home/agkanogiannis/software/vcftools_zlib_1.2.11/src/cpp/vcftools ;
bcftools=/home/agkanogiannis/bin/bcftools ;
vcfquery=/work/opt/vcftools/bin/vcf-query ;
plink=/home/agkanogiannis/software/plink2 ;
plink1=/home/agkanogiannis/software/plink ;
java=`which java` ;
pdftk=/home/agkanogiannis/software/pdftk/pdftk.jar ;
Rscript=/home/agkanogiannis/bin/Rscript ;

if [[ ! ${input_file} ]]; then
	usage ;
fi ;


temp_input=$(readlink -f ${input_file}) ;
mkdir -p ${basename} ;
pushd ${basename} ;
ln -sfn ${temp_input} . ;


####################################
## QC

### Step 1 ### Missingness
miss_in=${basename} ;

# Investigate missingness per individual and per SNP and make histograms
${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
	--vcf ${input_file} --missing --out ${miss_in} --silent ;
miss_in_v=$(echo "`wc -l ${miss_in}.vmiss|cut -d' ' -f1`-1"|bc -l) ;
miss_in_s=$(echo "`wc -l ${miss_in}.smiss|cut -d' ' -f1`-1"|bc -l) ;
${Rscript} --no-save <(echo -e "
	if (!require(reshape2)) install.packages('reshape2', repos='http://cran.us.r-project.org')
	if (!require(ggplot2)) install.packages('ggplot2', repos='http://cran.us.r-project.org')
	if (!require(lattice)) install.packages('lattice', repos='http://cran.us.r-project.org')
	suppressPackageStartupMessages(library(reshape2))
	suppressPackageStartupMessages(library(ggplot2))
	suppressPackageStartupMessages(library(lattice))
	args = commandArgs(trailingOnly=TRUE)
	snpmiss<-read.table(file=paste0(args[1],'.vmiss'),header=T,stringsAsFactors=F,fill=T,colClasses=c('NULL',NA,'NULL','NULL',NA))
	colnames(snpmiss) <- c('Variant','Variant Missingness')
	indmiss<-read.table(file=paste0(args[1],'.smiss'), header=T,stringsAsFactors=F,fill=T,colClasses=c(NA,'NULL','NULL',NA))
	colnames(indmiss) <- c('Sample','Sample Missingness')
	#pdf(file=paste0(paste0(args[1],'.miss.pdf')),paper='usr',width=11, height=8.5)
	vmiss.plot <-
		ggplot(snpmiss, aes(x=\`Variant Missingness\`)) +
	    geom_histogram(binwidth=0.01,fill='darkblue') +
	    labs(title='Variant Missingness \n(${miss_in_v} Variants)') +
	    scale_y_continuous(expand=c(0, 0)) +
	    scale_x_continuous(expand=c(0, 0), limits=c(0,1)) +
	    theme(plot.title=element_text(size=12))
	smiss.plot <-
		ggplot(indmiss, aes(x=\`Sample Missingness\`)) +
	    geom_histogram(binwidth=0.01,fill='darkblue') +
	    labs(title='Sample Missingness \n(${miss_in_s} Samples)') +
	    scale_y_continuous(expand=c(0, 0)) +
	    scale_x_continuous(expand=c(0, 0), limits=c(0,1)) +
	    theme(plot.title=element_text(size=12))
	save(vmiss.plot,smiss.plot,file='${basename}.1_miss.RData')
  	#print(vmiss.plot)
  	#print(smiss.plot)
  	#dev.off()
	") \
	${miss_in} ;

# Delete SNPs and individuals with high levels of missingness
if (( $(echo "${vmiss} < 1" | bc -l) && $(echo "${smiss} < 1" | bc -l) )); then 
	${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
		--vcf ${input_file} --geno ${vmiss} --export vcf-4.2 --out ${miss_in}.vmiss_${vmiss} --silent ;
	${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
		--vcf ${miss_in}.vmiss_${vmiss}.vcf --mind ${smiss} --export vcf-4.2 --out ${miss_in}.vmiss_${vmiss}.smiss_${smiss} --silent ;
	miss_out=${miss_in}.vmiss_${vmiss}.smiss_${smiss} ;
fi ;
if (( $(echo "${vmiss} < 1" | bc -l) && $(echo "${smiss} >= 1" | bc -l) )); then
	${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
		--vcf ${input_file} --geno ${vmiss} --export vcf-4.2 --out ${miss_in}.vmiss_${vmiss} --silent ;
	miss_out=${miss_in}.vmiss_${vmiss} ;
fi ;
if (( $(echo "${vmiss} >= 1" | bc -l) && $(echo "${smiss} < 1" | bc -l) )); then
	${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
		--vcf ${input_file} --mind ${smiss} --export vcf-4.2 --out ${miss_in}.smiss_${smiss} --silent ;
	miss_out=${miss_in}.smiss_${smiss} ;
fi ;
if (( $(echo "${vmiss} >= 1" | bc -l) && $(echo "${smiss} >= 1" | bc -l) )); then
	${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
		--vcf ${input_file} --export vcf-4.2 --out ${miss_in}.nofilt --silent ;
	miss_out=${miss_in}.nofilt ;
fi ;


### Step 2 ### MAF
maf_in=${miss_out} ;

# Generate a plot of the MAF distribution.
${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
	--vcf ${maf_in}.vcf --freq --out ${maf_in} --silent ;
maf_in_v=$(echo "`wc -l ${maf_in}.afreq|cut -d' ' -f1`-1"|bc -l) ;
maf_in_s=$(${vcfquery} -l ${maf_in}.vcf |wc -l|cut -d' ' -f1) ;
${Rscript} --no-save <(echo -e "
	suppressPackageStartupMessages(library(reshape2))
	suppressPackageStartupMessages(library(ggplot2))
	suppressPackageStartupMessages(library(lattice))
	args = commandArgs(trailingOnly=TRUE)
	afreq<-read.table(file=paste0(args[1],'.afreq'), header=T,stringsAsFactors=F,fill=T,colClasses=c('NULL',NA,'NULL','NULL',NA,'NULL'))
	colnames(afreq) <- c('Variant','MAF')
	afreq\$MAF <- as.numeric(afreq\$MAF)
	afreq <- afreq[!is.na(afreq\$MAF) & afreq\$MAF>0,]
	afreq\$MAF[afreq\$MAF>0.5] <- 1.0 - afreq\$MAF[afreq\$MAF>0.5]
	#pdf(file=paste0(paste0(args[1],'.afreq.pdf')),paper='usr',width=11, height=8.5)
	maf.plot <-
		ggplot(afreq, aes(x=MAF)) +
		labs(title='MAF \non V_miss=${vmiss} S_miss=${smiss} \n(${maf_in_v} Variants, ${maf_in_s} Samples)') +
	    geom_histogram(binwidth=0.01,fill='darkblue') +
	    scale_y_continuous(expand=c(0, 0)) +
	    theme(panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),axis.text.x=element_text(angle=45),plot.title=element_text(size=12)) + 
	    scale_x_log10(breaks=c(0.001,0.002,0.003,0.004,0.005,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,1.0),limits=NULL)
	save(maf.plot,file='${basename}.2_maf.RData')
  	#print(maf.plot)
  	#dev.off()
	") \
	${maf_in} ;

# Remove SNPs with a low and high MAF frequency.
if (( $(echo "${maf} > 0" | bc -l) )); then
	${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
		--vcf ${maf_in}.vcf --maf ${maf} minor --max-maf ${maf_max} minor --export vcf-4.2 --out ${maf_in}.MAF_${maf}_${maf_max} --silent ;
	maf_out=${maf_in}.MAF_${maf}_${maf_max} ;
else
	maf_out=${maf_in} ;
fi ;


### Step 3 ### HWE
hwe_in=${maf_out} ;

# Check the distribution of HWE Fisher exact test p-values of all SNPs.
${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
	--vcf ${hwe_in}.vcf --hardy --out ${hwe_in} --silent ;
hwe_in_v=$(echo "`wc -l ${hwe_in}.hardy|cut -d' ' -f1`-1"|bc -l) ;
hwe_in_s=$(${vcfquery} -l ${hwe_in}.vcf |wc -l|cut -d' ' -f1) ;
${Rscript} --no-save <(echo -e "
	suppressPackageStartupMessages(library(reshape2))
	suppressPackageStartupMessages(library(ggplot2))
	suppressPackageStartupMessages(library(lattice))
	args = commandArgs(trailingOnly=TRUE)
	hwe<-read.table(file=paste0(args[1],'.hardy'), header=T,stringsAsFactors=F,fill=T,colClasses=c('NULL',NA,'NULL','NULL','NULL','NULL','NULL','NULL','NULL',NA))
	colnames(hwe) <- c('Variant','HWE p-value')
	hwe\$\`HWE p-value\` <- as.numeric(as.character(hwe\$\`HWE p-value\`))
	hwe <- hwe[!is.na(hwe\$\`HWE p-value\`) & hwe\$\`HWE p-value\`>0,]
	#pdf(file=paste0(paste0(args[1],'.hardy.pdf')),paper='usr',width=11, height=8.5)
	hwe.plot <-
		ggplot(hwe, aes(x=\`HWE p-value\`)) +
		labs(title='HWE p-value \non V_miss=${vmiss} S_miss=${smiss} MAF=${maf} \n(${hwe_in_v} Variants, ${hwe_in_s} Samples)') +
	    geom_histogram(binwidth=0.01,fill='darkblue') +
	    scale_y_continuous(expand=c(0, 0)) +
	    theme(panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),plot.title=element_text(size=12)) +
	    scale_x_continuous(expand=c(0, 0), limits=c(0,1))
	hwe_zoom.plot <-
		ggplot(subset(hwe,\`HWE p-value\`<1e-3), aes(x=\`HWE p-value\`)) +
		labs(title='HWE p-value (low end)\non V_miss=${vmiss} S_miss=${smiss} MAF=${maf} \n(${hwe_in_v} Variants, ${hwe_in_s} Samples)') +
	    geom_histogram(binwidth=0.05,fill='darkblue') +
	    scale_y_continuous(expand=c(0, 0)) +
	    theme(panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),axis.text.x=element_text(angle=45,size=5),plot.title=element_text(size=12)) +
	    scale_x_log10(breaks=c(1e-30,1e-20,1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3))
	save(hwe.plot,hwe_zoom.plot,file='${basename}.3_hwe.RData')
  	#print(hwe.plot)
  	#print(hwe_zoom.plot)
  	#dev.off()
	") \
	${hwe_in} ;

# The HWE threshold filters out only SNPs which deviate extremely from HWE.
if (( $(echo "${hwe_tmp} > 0" | bc -l) )); then
	${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
		--vcf ${hwe_in}.vcf --hwe ${hwe} --export vcf-4.2 --out ${hwe_in}.hardy_${hwe} --silent ;
	hwe_out=${hwe_in}.hardy_${hwe} ;
else
	hwe_out=${hwe_in} ;
fi ;


### Step 4 ### LD pruned snps
ld_in=${hwe_out} ;

# To generate a list of non-(highly)correlated SNPs and prune the SNPs using the command --indep-pairwise’.
# The parameters ‘50 5 0.2’ stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.
if (( `echo "${ld}>=0 && ${ld}<1"|bc -l` )) ; then
	${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
		--vcf ${ld_in}.vcf --indep-pairwise 50 5 ${ld} --out ${ld_in} --silent ;
	# Save LD snps
	${plink} --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
		--vcf ${ld_in}.vcf --extract ${ld_in}.prune.in --export vcf-4.2 --out ${ld_in}.ld --silent ;
	ld_out=${ld_in}.ld ;
else
	ld_out=${ld_in} ;
fi ;


### Step 5 ### Heterozygosity (on LD pruned snps)
het_in=${ld_out} ;

# Generate a plot of the distribution of the heterozygosity rate of your subjects.
# And save a list of individuals with a heterozygosity rate deviating more than 4 sd from the mean.

${plink1} --keep-allele-order --real-ref-alleles --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
	--vcf ${het_in}.vcf --het --out ${het_in} --silent ;
het_in_v=$(${bcftools} stats ${het_in}.vcf|grep 'number of records:'|cut -f4) ;
het_in_s=$(${vcfquery} -l ${het_in}.vcf |wc -l|cut -d' ' -f1) ;

# Plot of the heterozygosity rate distribution
${Rscript} --no-save <(echo -e "
	suppressPackageStartupMessages(library(reshape2))
	suppressPackageStartupMessages(library(ggplot2))
	suppressPackageStartupMessages(library(lattice))
	args = commandArgs(trailingOnly=TRUE)
	het<-read.table(file=paste0(args[1],'.het'), header=T,stringsAsFactors=F,fill=T)
	het\$\`Heterozygosity\` <- (het\$\`N.NM.\` - het\$\`O.HOM.\`)/het\$\`N.NM.\`
	het\$\`Heterozygosity\` <- as.numeric(as.character(het\$\`Heterozygosity\`))
	#pdf(file=paste0(paste0(args[1],'.het.pdf')),paper='usr',width=11, height=8.5)
	het.plot <-
		ggplot(het, aes(x=\`Heterozygosity\`)) +
		labs(title='Sample Heterozygosity \non V_miss=${vmiss} S_miss=${smiss} MAF=${maf} HWE=${hwe} \n(LD=${ld}, ${het_in_v} Variants, ${het_in_s} Samples)') +
	    geom_histogram(binwidth=0.01,fill='darkblue')+
	    scale_y_continuous(expand=c(0, 0)) +
	    scale_x_continuous(expand=c(0, 0), limits=c(0,1)) +
	    theme(plot.title=element_text(size=12))
	f.plot <-
		ggplot(het, aes(x=F)) +
		labs(title='Sample Inbreeding \non V_miss=${vmiss} S_miss=${smiss} MAF=${maf} HWE=${hwe} \n(LD=${ld}, ${het_in_v} Variants, ${het_in_s} Samples)') +
	    geom_histogram(binwidth=0.01,fill='darkblue')+
	    scale_y_continuous(expand=c(0, 0)) +
	    scale_x_continuous(expand=c(0, 0), limits=c(-1,1)) +
	    theme(plot.title=element_text(size=12))
	save(het.plot,f.plot,file='${basename}.4_het.RData')
  	#print(het.plot)
  	#print(f.plot)
  	#dev.off()
	") \
	${het_in} ;

# The following code generates a list of individuals who deviate more than 4 standard deviations from the heterozygosity rate mean.
${Rscript} --no-save <(echo -e "
	args = commandArgs(trailingOnly=TRUE)
	het <- read.table(paste0(args[1],'.het'), head=TRUE)
	het\$HET_RATE <- (het\$'N.NM.' - het\$'O.HOM.')/het\$'N.NM.'
	het_fail <- subset(het, (het\$HET_RATE < mean(het\$HET_RATE)-4*sd(het\$HET_RATE)) | (het\$HET_RATE > mean(het\$HET_RATE)+4*sd(het\$HET_RATE)))
	het_fail\$HET_DST <- (het_fail\$HET_RATE-mean(het\$HET_RATE))/sd(het\$HET_RATE)
	write.table(het_fail, paste0(args[1],'.het_fail'), row.names=FALSE)") \
	${het_in} ;
# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the IID column
sed 's/"// g' ${het_in}.het_fail | awk '{print $2}'> ${het_in}.het_fail_ind ;
sed -i '1d' ${het_in}.het_fail_ind ;
sed -i '1i #IID' ${het_in}.het_fail_ind ;

het_out=${het_in} ;

####################################

# vcf INFO stats
####################################
if [[ "${info_stats}" -eq 1 ]]; then
	${Rscript} --no-save <(echo -e "
		if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='http://cran.us.r-project.org')
		if (!require(VariantAnnotation)) BiocManager::install('VariantAnnotation')
		suppressPackageStartupMessages(library(VariantAnnotation))
		suppressPackageStartupMessages(library(ggplot2))
		args = commandArgs(trailingOnly=TRUE)
		vcf<-args[1]
		INFO <- data.frame(QD=readInfo(file=vcf,x='QD'),
                   		   MQ=readInfo(file=vcf,x='MQ'),
                   		   FS=readInfo(file=vcf,x='FS'),
                   		   SOR=readInfo(file=vcf,x='SOR'),
                   		   MQRankSum=readInfo(file=vcf,x='MQRankSum'),
                   		   ReadPosRankSum=readInfo(file=vcf,x='ReadPosRankSum'))
        QD_p <- ggplot(INFO, aes(x=\`QD\`)) + geom_density(color='darkblue', fill='lightblue')
        MQ_p <- ggplot(INFO, aes(x=\`MQ\`)) + geom_density(color='darkblue', fill='lightblue')
		FS_p <- ggplot(INFO, aes(x=\`FS\`)) + geom_density(color='darkblue', fill='lightblue') + scale_x_log10(limits=NULL)
		SOR_p <- ggplot(INFO, aes(x=\`SOR\`)) + geom_density(color='darkblue', fill='lightblue')
		MQR_p <- ggplot(INFO, aes(x=\`MQRankSum\`)) + geom_density(color='darkblue', fill='lightblue')
		RPR_p <- ggplot(INFO, aes(x=\`ReadPosRankSum\`)) + geom_density(color='darkblue', fill='lightblue')
        save(QD_p,MQ_p,FS_p,SOR_p,MQR_p,RPR_p,file='${basename}.5_INFO_stats.RData')
		") \
		${ld_out}.vcf ;
fi ;
####################################

####################################
## Relatedness (on LD pruned snps)

rel_in=${ld_out} ;

# Assuming a random population sample we are going to exclude all individuals above the pihat threshold
# Save clone and related lists *.clones.tsv *.close_related.tsv *.clusters.tsv
# ld=1.0 as we use the already ld pruned snp set
${Rscript} --vanilla ${relatedness_R} ${rel_in}.vcf ${pihat} 1.0 ${MLE} ${pop_groups} ;

rel_out=${rel_in} ;

####################################


# Collate QC graphs
${Rscript} --no-save <(echo -e "
	if (!require(gridExtra)) install.packages('gridExtra', repos='http://cran.us.r-project.org')
	suppressPackageStartupMessages(library(gridExtra))
	pdf(file='${basename}.QC.report.pdf',paper='usr',width=11, height=8.5)
	load(file='${basename}.1_miss.RData')
	load(file='${basename}.2_maf.RData')
	load(file='${basename}.3_hwe.RData')
	load(file='${basename}.4_het.RData')
	grid.arrange(vmiss.plot,smiss.plot,maf.plot,hwe.plot,hwe_zoom.plot,het.plot,f.plot,
 		widths = c(1,1,2),
  		layout_matrix=rbind(c(1,2,3),c(4,4,5),c(6,6,7)))
	dev.off()
	file.remove('${basename}.1_miss.RData')
	file.remove('${basename}.2_maf.RData')
	file.remove('${basename}.3_hwe.RData')
	file.remove('${basename}.4_het.RData')
	") \
	;


# Collate INFO_stats graphs
if [[ "${info_stats}" -eq 1 ]]; then
	${Rscript} --no-save <(echo -e "
		if (!require(gridExtra)) install.packages('gridExtra', repos='http://cran.us.r-project.org')
		suppressPackageStartupMessages(library(gridExtra))
		pdf(file='${basename}.INFO.stats.pdf',paper='usr',width=11, height=8.5)
		load(file='${basename}.5_INFO_stats.RData')
		grid.arrange(QD_p,MQ_p,FS_p,SOR_p,MQR_p,RPR_p,
	 		widths = c(1,1),
	  		layout_matrix=rbind(c(1,2),c(3,4),c(5,6)))
		dev.off()
		file.remove('${basename}.5_INFO_stats.RData')
		") \
		;
fi ;


# Collate QC, INFO_stats and relatedness graphs

if [[ "${info_stats}" -eq 1 ]]; then
	${java} -jar ${pdftk} ${basename}.QC.report.pdf ${basename}.INFO.stats.pdf ${rel_in}.vcf.relatedness.pdf cat output ${basename}.QC_relatedness.report.pdf ;
	rm -f  ${basename}.QC.report.pdf ${basename}.INFO.stats.pdf ${rel_in}.vcf.relatedness.pdf Rplots.pdf ;
else
	${java} -jar ${pdftk} ${basename}.QC.report.pdf ${rel_in}.vcf.relatedness.pdf cat output ${basename}.QC_relatedness.report.pdf ;
	rm -f  ${basename}.QC.report.pdf ${rel_in}.vcf.relatedness.pdf Rplots.pdf ;
fi ;

ln -sfn ${hwe_out}.vcf 				${basename}.initialFiltered.vcf ;
ln -sfn  ${ld_out}.vcf 				${basename}.initialFiltered.ld.vcf ;
ln -sfn ${rel_out}.vcf.ibd_mom.dist ${basename}.initialFiltered.vcf.ibd_mom.dist ;

popd ;

exit 0;
