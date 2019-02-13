enigma_tsv=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/data/known_variants/enigma_variants_20181221.tsv
enigma_txt=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/data/known_variants/enigma_variants_20181221.txt

## Use this to get col numbers
#head -1 enigma_variants_20181221.tsv | sed 's/\t/\n/g' | nl

# 1	id
# 15	Clinical_significance_ENIGMA
# 111	Hg37_Start ## Use this as Pos instead
# 112	Hg37_End
# 122	Chr
# 123	Pos ## Pos is based on hg19 (Hg_38 in tsv)
# 124	Ref
# 125	Alt

# echo -e "id\tchrom\tpos\tref\talt\tsig" > $enigma_txt
# cat $enigma_tsv | grep -vE '^#' | awk '{
# 	chrom=$122
# 	pos=$111
# 	id=$1
# 	ref=$124
# 	alt=$125
# 	sig=gsub(" ","_", $15)
	
# 	print id "\t" chrom "\t" pos "\t" ref "\t" alt "\t" sig
	
# }' >> $enigma_txt

## awk doesn't work. Use excel instead
