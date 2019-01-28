clinvar_vcf=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/data/variant_significance/clinvar_20181217.vcf.gz
clinvar_txt=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/data/variant_significance/clinvar_20181217.txt

echo -e "id\tchrom\tpos\tref\talt\tsig" > $clinvar_txt
zcat $clinvar_vcf | grep -vE '^#' | while read chrom pos id ref alt qual filter info; do
	sig=$(echo $info | grep -oP "CLNSIG=[a-zA-z_]+" | sed 's/CLNSIG=//g')
	if [[ -z $sig ]]; then sig="blank"; fi
	echo -e "${id}\t${chrom}\t${pos}\t${ref}\t${alt}\t${sig}"
done > $clinvar_txt

## Deprecated
# echo -e "id\tchrom\tpos\tref\talt\tsig" > $clinvar_txt
# zcat $clinvar_vcf | grep -vE '^#' | awk '{
# 	chrom=$1
# 	pos=$2
# 	id=$3
# 	ref=$4
# 	alt=$5

# 	match($0, "CLNSIG=[a-zA-z_]+", sig)
# 	gsub("CLNSIG=","", sig[0])
	
# 	print id "\t" chrom "\t" pos "\t" ref "\t" alt "\t" sig[0]
	
# }' >> $clinvar_txt