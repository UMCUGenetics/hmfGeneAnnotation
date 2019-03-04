#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/
source $ROOT_DIR/loadPaths.sh

subsetCadd (){
	tsv_in=$1
	#tsv_in=$ROOT_DIR/data/variant_significance/cadd/whole_genome_SNVs_inclAnno.tsv.gz
	
	tsv_ss=$(dirname $tsv_in)/$(basename $tsv_in '.tsv.gz')_ss.tsv.gz
	#tsv_ss2=$(dirname $tsv_in)/$(basename $tsv_in '.tsv.gz')_ss2.tsv.gz

	## Subset
	zcat $tsv_in | $JAVA -jar $SNPSIFT intervals $GENES_BED | ## Extract only variants that are in selected genes
	grep -vE '^##' | ## Remove vcf header lines added by SnpSift
	awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$(NF-1),"\t",$NF}' | uniq -u | ## Extract relevant cols, and unique rows
	gzip -c > $tsv_ss

}

# if [[ ! -f $tsv_ss2 ]]; then
# 	zcat $tsv_ss1 |
# 	awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$(NF-1),"\t",$NF}' | 
# 	uniq -u | gzip -c > $tsv_ss2
# fi

## Extract only variants that are in selected genes
# if [[ ! -f $tsv_ss2 ]]; then
# 	zcat $tsv_ss1 | grep '#Chr' | head -n 1 | gzip -c > $tsv_ss2

# 	bed_chrom=($(cat $GENES_BED | grep -vE '^#' | awk '{print $1}'))
# 	bed_start=($(cat $GENES_BED | grep -vE '^#' | awk '{print $2}'))
# 	bed_end=($(cat $GENES_BED | grep -vE '^#' | awk '{print $3}'))

# 	zcat $tsv_ss1 | grep -vE '^#' | while read chrom pos etc; do
# 		for ((i=0; i<${#bed_chrom[@]}; ++i)); do
# 			if [[ ${bed_chrom[i]} -eq $chrom && ${bed_start[i]} -le $pos && ${bed_end[i]} -ge $pos ]]; then
# 				echo -e "${chrom}\t${pos}\t${etc}"
# 				break
# 			fi
# 		done
# 	done | gzip -c >> $tsv_ss2
# fi

#zcat $tsv_in | $JAVA -jar $SNPSIFT intervals $GENES_BED | grep -vE '^##' | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$(NF-1),"\t",$NF}' | uniq -u


# tsv_ss1=$ROOT_DIR/data/variant_significance/cadd/InDels_inclAnno_ss1.tsv.gz
# #tsv_ss1=$ROOT_DIR/data/variant_significance/cadd/whole_genome_SNVs_inclAnno_ss1.tsv.gz

# bed_chrom=($(cat $GENES_BED | grep -vE '^#' | awk '{print $1}'))
# bed_start=($(cat $GENES_BED | grep -vE '^#' | awk '{print $2}'))
# bed_end=($(cat $GENES_BED | grep -vE '^#' | awk '{print $3}'))

# # x=1314891
# # for ((i=0; i<${#bed_chrom[@]}; ++i)); do
# # 	if [[ ${bed_start[i]} -le $x && ${bed_end[i]} -ge $x ]]; then echo ${bed_start[i]}; fi
# # done

# zcat $tsv_ss1 | grep -vE '^#' | while read chrom pos etc; do
	
# 	for ((i=0; i<${#bed_chrom[@]}; ++i)); do
# 		if [[ ${bed_chrom[i]} -eq $chrom && ${bed_start[i]} -le $pos && ${bed_end[i]} -ge $pos ]]; then
# 			echo -e "${chrom}\t${pos}\t${etc}"
# 			break
# 		fi
# 	done
# done

#zcat $tsv_ss1 | $JAVA -jar $SNPSIFT intervals $GENES_BED | grep -vE '^##' | gzip -c > $tsv_ss2

