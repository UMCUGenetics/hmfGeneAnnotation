#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

source $ROOT_DIR/loadPaths.sh

#========= Filter vcfs  =========#
## by PASS and selected gene coords. For germline vcf, also remove somatic variants
filterVcf(){
	vcf_in=$1
	vcf_out=$2
	genes_bed=$3
	mode=$4

	if [[ $mode == "somatic" ]]; then
		filt_string="(FILTER='PASS')"
	elif [[ $mode == "germline" ]]; then
		filt_string="(FILTER='PASS') & !(GEN[0].GT='0/0') & !(GEN[0].GT='./.')"
	else 
		echo Available modes are: \"somatic\", \"germline\"
	fi

	zcat $vcf_in |
	$JAVA -jar $SNPSIFT intervals $genes_bed |
	$JAVA -jar $SNPSIFT filter "$filt_string" |
	gzip -c > $vcf_out
}
