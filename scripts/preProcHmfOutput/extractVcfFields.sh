#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

source $ROOT_DIR/loadPaths.sh

extractVcfFields(){
	vcf_in=$1
	txt_out=$2
	mode=$3

	# vcf_in=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020719R_CPCT02020719T/CPCT02020719R_CPCT02020719T.germ.vcf.gz
	# txt_out=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020719R_CPCT02020719T/CPCT02020719R_CPCT02020719T.germ.txt.gz
	# mode='germline'

	if [[ $mode == "somatic" ]]; then
		tumor_col=0
	elif [[ $mode == "germline" ]]; then
		tumor_col=1
	else 
		echo Available modes are: \"somatic\", \"germline\"
	fi

	$JAVA -jar $SNPSIFT extractFields $vcf_in \
		CHROM POS REF ALT \
		ANN[0].EFFECT \
		ANN[0].GENE \
		ANN[0].GENEID \
		ANN[0].FEATUREID \
		ANN[0].HGVS_C \
		GEN[$tumor_col].GT \
		GEN[$tumor_col].AD[0] GEN[$tumor_col].AD[1] | \
	sed -e "1s/.*/chrom\tpos\tref\talt\tsnpeff_eff\tsnpeff_gene\tensembl_gene_id\tensembl_transcript_id\thgvs_c\ttumor_gt\ttumor_ad_ref\ttumor_ad_alt/" | \
	gzip -c > $txt_out
}

# extractVcfFields \
# /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020459R_CPCT02020459T/CPCT02020459R_CPCT02020459T.germ.vcf.gz \
# /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020459R_CPCT02020459T/CPCT02020459R_CPCT02020459T.germ.txt.gz \
# 'germline'