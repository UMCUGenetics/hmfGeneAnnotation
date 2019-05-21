#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

source $ROOT_DIR/loadPaths.sh

extractVcfFields(){
	vcf_in=$1
	txt_out=$2
	mode=$3

	## Snpsift
	if [[ $mode == 'ignore' ]]; then
		GEN_GT=''
		GEN_AD0=''
		GEN_AD1=''

	elif [[ $mode == 'somatic' ]]; then
		GEN_GT='GEN[0].GT'
		GEN_AD0='GEN[0].AD[0]'
		GEN_AD1='GEN[0].AD[1]'

	elif [[ $mode == 'germline' ]]; then
		GEN_GT='GEN[1].GT'
		GEN_AD0='GEN[1].AD[0]'
		GEN_AD1='GEN[1].AD[1]'
	else 
		echo Available modes are: \"somatic\", \"germline\", \"ignore\"
	fi

	$JAVA -jar $SNPSIFT extractFields $vcf_in \
		CHROM POS REF ALT \
		ANN[0].EFFECT \
		ANN[0].GENE \
		ANN[0].GENEID \
		ANN[0].FEATUREID \
		ANN[0].HGVS_C \
		$GEN_GT \
		$GEN_AD0 \
		$GEN_AD1 |
		gzip -c > ${txt_out}.temp

	
	## Insert NAs for columns with no values
	## Create header
	zcat ${txt_out}.temp | 
	
	if [[ $mode == 'ignore' ]]; then
		awk '{print $0"\t""NA""\t""NA""\t""NA"}'
	else
		awk '{print $0}'
	fi |
	
	sed -e "1s/.*/chrom\tpos\tref\talt\tsnpeff_eff\tsnpeff_gene\tensembl_gene_id\tensembl_transcript_id\thgvs_c\ttumor_gt\ttumor_ad_ref\ttumor_ad_alt/" |
	gzip -c > $txt_out && rm ${txt_out}.temp

	
}

# extractVcfFields \
# /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020459R_CPCT02020459T/CPCT02020459R_CPCT02020459T.germ.vcf.gz \
# /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020459R_CPCT02020459T/CPCT02020459R_CPCT02020459T.germ.txt.gz \
# 'germline'