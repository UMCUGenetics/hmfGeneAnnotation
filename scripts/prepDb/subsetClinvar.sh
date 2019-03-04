#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/
source $ROOT_DIR/loadPaths.sh

clinvar_vcf=$ROOT_DIR/data/variant_significance/clinvar/clinvar_20181217.vcf.gz

## Output files
clinvar_ss_vcf=$ROOT_DIR/data/variant_significance/clinvar/clinvar_20181217_ss.vcf.gz
clinvar_ss_txt=$ROOT_DIR/data/variant_significance/clinvar/clinvar_20181217_ss.txt

## Subset variants by provided gene list to optimize for speed
zcat $clinvar_vcf | $JAVA -jar $SNPSIFT intervals $GENES_BED | gzip -c > $clinvar_ss_vcf

## Extract relevant fields
echo -e "chrom\tpos\tref\talt\tsig\tid" > $clinvar_ss_txt
zcat $clinvar_ss_vcf | grep -vE '^#' | while read chrom pos id ref alt qual filter info; do
	sig=$(echo $info | grep -oP "CLNSIG=[a-zA-z_]+" | sed 's/CLNSIG=//g')
	if [[ -z $sig ]]; then sig="blank"; fi
	echo -e "${chrom}\t${pos}\t${ref}\t${alt}\t${sig}\t${id}"
done > $clinvar_ss_txt && rm $clinvar_ss_vcf
