#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/
source $ROOT_DIR/loadPaths.sh

detIsHotspotMut (){
	in_txt=$1
	out_txt=$2

	echo 'is_hotspot_mut' | gzip -c > $out_txt
	zcat $in_txt | tail -n +2 | while read chrom pos ref alt etc; do

		## Debugging
		#has_match=$(zcat $HOTSPOTS_DB | grep -cP -m1 "^${chrom}\t${pos}\t${ref}\t${alt}")
		#echo -e "${has_match}\t${chrom}\t${pos}\t${ref}\t${alt}"
		
		zgrep -cP -m1 "^${chrom}\t${pos}\t${ref}\t${alt}" $HOTSPOTS_DB
	done | gzip -c >> $out_txt
}


## Test
# wd=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020493R_CPCT02020493T/
# in_txt=$wd/CPCT02020493R_CPCT02020493T.germ.txt.gz
# out_txt=$wd/varsig/hotspots_germ.txt.gz

# wd=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/data/variant_significance/HMF_hotspots/
# in_txt=$wd/CPCT02020493R_CPCT02020493T.germ.txt.gz
# out_txt=$wd/hotspots_germ.txt.gz

# detIsHotspotMut $in_txt $out_txt