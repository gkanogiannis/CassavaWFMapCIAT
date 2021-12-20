#!/bin/bash

#
# README_geno_preprocess.txt.sh
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


mother=ECU72 ;
father=TMS60444 ;
base_work_folder=/data_A/Cassava/Cassava_WGS/Cassava_WGS_CIAT_GM8586_MAP/vcf_work ;
base_name=GM8586 ;
p_value=0.05 ;


# * 1.select_SNPs
mkdir -p ${base_work_folder}/1.select_SNPs ;
pushd ${base_work_folder}/1.select_SNPs ;

# Remove erroneous and useless SNPs (inconsistent child segregation, both parents hom, missing parents)
vcf_initial=${base_name}.initial ;
groovy /home/agkanogiannis/CIAT/scripts/FilterVCFforMap.groovy <(zcat /data_A/Cassava/Cassava_WGS/Cassava_WGS_CIAT_GM8586_MAP/vcf_unique/vcf/Cassava_WGS_CIAT_GM8586_MAP.snps.filter_info_minDP10_maxDP40_GQ30.vcf.gz) ${mother} ${father} ${p_value} 2>/dev/null > ${vcf_initial}.vcf

# Run first pass QC and filtering
#bash /home/agkanogiannis/CIAT/scripts/initialFilteringAndRelatedness.sh -i ${vcf_initial}.vcf ;

# Investigate first pass and run second pass
bash /home/agkanogiannis/CIAT/scripts/initialFilteringAndRelatedness.sh -i ${vcf_initial}.vcf -f -v 0.10 -s 0.20 ;
pushd ${vcf_initial} ;
vcf_filtered=${vcf_initial}.initialFiltered ;

# Plot distances between crosses and parents
Rscript /home/agkanogiannis/CIAT/scripts/plot_dist_parents.R ${vcf_filtered}.vcf.ibd_mom.dist ${father} ${mother} ;

# Calculate vcftools relatedness
# Plot relatedness between crosses and parents
vcftools --vcf ${vcf_filtered}.vcf --relatedness2 --out ${vcf_filtered} ;
rm ${vcf_filtered}.log ;
Rscript /home/agkanogiannis/CIAT/scripts/plot_relatedness_parents.R ${vcf_filtered}.relatedness2 ${father} ${mother} ;

# Remove problematic 
# Remove clones
# GM8586-50, GM8586-36, GM8586-72
popd ;
vcf_final=${base_name}.final ;
vcftools \
	--vcf ${vcf_initial}/${vcf_filtered}.vcf \
	--remove-indv GM8586-50 --remove-indv GM8586-36 --remove-indv GM8586-72 \
	--recode --recode-INFO-all \
	--stdout > ${vcf_final}.vcf ;
	# --exclude-positions <(cat ${vcf_initial}/${vcf_nosegr}.positions) \
	# --exclude-positions <(cat ${vcf_initial}/${vcf_crosses_nomend}.positions_vcftools ${vcf_initial}/${vcf_crosses_nosegr}.positions ${vcf_initial}/${vcf_nosegr}.positions) \

# Plot distances between crosses and parents
bash /home/agkanogiannis/CIAT/scripts/initialFilteringAndRelatedness.sh -f -i ${vcf_final}.vcf ;
Rscript /home/agkanogiannis/CIAT/scripts/plot_dist_parents.R ${vcf_final}/${vcf_final}.nofilt.vcf.ibd_mom.dist ${father} ${mother} ;

# Calculate relatedness
# Plot relatedness between crosses and parents
pushd ${vcf_final} ;
vcftools --vcf ${vcf_final}.vcf --relatedness2 --out ${vcf_final} ;
Rscript /home/agkanogiannis/CIAT/scripts/plot_relatedness_parents.R ${vcf_final}.relatedness2 ${father} ${mother} ;
popd ;


#* 2.convert_SNPs_JoinMap
popd ;
mkdir -p ${base_work_folder}/2.convert_SNPs_JoinMap ;
pushd ${base_work_folder}/2.convert_SNPs_JoinMap ;

vcf_joinmap=GM8586 ;
ln -sfn ../encoding-snp.py . ;
ln -sfn ../1.select_SNPs/${vcf_final}.vcf ${vcf_joinmap}.vcf ;

# Convert vcf to joinmap format (use Vianey script)
python2.7 \
	encoding-snp.py \
	-i ${vcf_joinmap}.vcf -m ${mother} -f ${father} \
	-o ${vcf_joinmap} ;

# Select <nnxnp>, <lmxll> and <hkxhk> SNPs
for i in nnxnp lmxll hkxhk; do
	vcftools \
		--vcf ${vcf_joinmap}.vcf \
		--recode --recode-INFO-all \
		--positions <(grep "$i" ${vcf_joinmap}.loc | awk '{print $1}' | sed 's/_/\t/g; s/^C/Chromosome/g; s/^S/Scaffold/g;') \
		--stdout > ${vcf_joinmap}.${i}.vcf ;
done ;

# LD pruning
for i in nnxnp lmxll hkxhk; do
	plink --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
                --vcf ${vcf_joinmap}.${i}.vcf --indep-pairwise 50 5 0.995 --out ${vcf_joinmap}.${i} 1>/dev/null ;
        plink --allow-extra-chr --set-missing-var-ids '@_#' --const-fid \
                --vcf ${vcf_joinmap}.${i}.vcf --extract ${vcf_joinmap}.${i}.prune.in --export vcf-4.2 --out ${vcf_joinmap}.${i}.ld 1>/dev/null ;
	bgzip ${vcf_joinmap}.${i}.ld.vcf ;
	tabix -p vcf ${vcf_joinmap}.${i}.ld.vcf.gz ;
done ;
vcf-concat ${vcf_joinmap}.nnxnp.ld.vcf.gz ${vcf_joinmap}.lmxll.ld.vcf.gz ${vcf_joinmap}.hkxhk.ld.vcf.gz > ${vcf_joinmap}.ld.vcf ;
gunzip *.gz ;
rm -rf *.tbi ; 

# Convert vcf to joinmap format
python2.7 \
	encoding-snp.py \
	-i ${vcf_joinmap}.ld.vcf -m ${mother} -f ${father} \
	-o ${vcf_joinmap}.ld ;
for i in nnxnp lmxll hkxhk; do
	python2.7 \
		encoding-snp.py \
		-i ${vcf_joinmap}.${i}.vcf -m ${mother} -f ${father} \
		-o ${vcf_joinmap}.${i} ;
	python2.7 \
		encoding-snp.py \
		-i ${vcf_joinmap}.${i}.ld.vcf -m ${mother} -f ${father} \
		-o ${vcf_joinmap}.${i}.ld ;
done ;
