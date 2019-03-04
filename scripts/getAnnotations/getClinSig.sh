#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/
source $ROOT_DIR/loadPaths.sh

getClinSig (){
	in_txt=$1
	out_txt=$2

	echo -e "clinvar_sig\tenigma_sig" | gzip -c > $out_txt
	zcat $in_txt | tail -n +2 | while read chrom pos ref alt etc; do
		
		## Retrieve pathogenicity
		clinvar_sig=$(grep -m 1 -P "^${chrom}\t${pos}\t${ref}\t${alt}" $CLINVAR_DB | cut -f 5 )
		enigma_sig=$(grep -m 1 -P "^${chrom}\t${pos}\t${ref}\t${alt}" $ENIGMA_DB | cut -f 5 )


		if [[ -z $clinvar_sig ]]; then clinvar_sig=NA; fi
		if [[ -z $enigma_sig ]]; then enigma_sig=NA; fi

		echo -e "${clinvar_sig}\t${enigma_sig}"
	done | gzip -c >> $out_txt
}

# ## Exec
# wd=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020386R_CPCT02020386T/
# in_txt=$wd/CPCT02020386R_CPCT02020386T.som.txt.gz
# out_txt=$wd/varsig_som/enigma.txt.gz
# db='enigma'

# getClinSig $in_txt $out_txt $db



#command=paste; for i in *.gz; do command="$command <(gzip -cd $i)"; done; eval $command
