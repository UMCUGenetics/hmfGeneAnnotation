#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

source $ROOT_DIR/loadPaths.sh

#========= Filter vcfs  =========#
## by PASS and selected gene coords. For germline vcf, also remove somatic variants
filterVcf(){
	vcf_in=$1
	vcf_out=$2
	mode=$3

	# vcf_in=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/180426_HMFregCPCT_FR13997274_FR16982076_CPCT02020719/180426_HMFregCPCT_FR13997274_FR16982076_CPCT02020719.annotated.vcf.gz
	# vcf_out=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020719R_CPCT02020719T/CPCT02020719R_CPCT02020719T.germ.vcf.gz
	# mode='germline'
	# bed=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/data/gene_selection/genes.bed

	if [[ $mode == "somatic" ]]; then
		filt_string="(FILTER='PASS')"
	elif [[ $mode == "germline" ]]; then
		filt_string="(FILTER='PASS') & !(GEN[0].GT='0/0') & !(GEN[0].GT='./.')"
	else 
		echo Available modes are: \"somatic\", \"germline\"
	fi

	zcat $vcf_in |
	$JAVA -jar $SNPSIFT intervals $GENES_BED |
	$JAVA -jar $SNPSIFT filter "$filt_string" |
	gzip -c > $vcf_out
}

# filterVcfs \
# /hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/170320_HMFregCPCT_FR13274499_FR10244814_CPCT02020459/170320_HMFregCPCT_FR13274499_FR10244814_CPCT02020459.annotated.vcf.gz \
# /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020459R_CPCT02020459T/CPCT02020459R_CPCT02020459T.germ.vcf.gz \
# 'germline'