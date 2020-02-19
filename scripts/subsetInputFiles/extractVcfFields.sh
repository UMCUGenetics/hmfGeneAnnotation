#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

source $ROOT_DIR/loadPaths.sh

extractVcfFields(){
	vcf_in=$1
	txt_out=$2

	$JAVA -jar $SNPSIFT extractFields $vcf_in \
		CHROM POS REF ALT \
		ANN[0].EFFECT \
		ANN[0].GENE \
		ANN[0].GENEID \
		ANN[0].FEATUREID \
		ANN[0].HGVS_C |
		gzip -c > ${txt_out}.temp

	## Edit header
	zcat ${txt_out}.temp |
	sed -e "1s/.*/chrom\tpos\tref\talt\tsnpeff_eff\tsnpeff_gene\tensembl_gene_id\tensembl_transcript_id\thgvs_c/" |
	gzip -c > $txt_out && rm ${txt_out}.temp
}