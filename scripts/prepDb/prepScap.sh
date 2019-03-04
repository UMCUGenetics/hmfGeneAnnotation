#!/bin/bash

module load tabix/1.7 ## tabix contains bgzip module


base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

scap_dir=$base_dir/data/variant_significance/SCAP/
raw_dir=$scap_dir/raw/

echo -e "\n#========= cbind splice type and merge to one file =========#"
addSpliceType (){
	in_file=$1
	type=$2
	
	zcat $in_file | tail -n +2 | awk -v type="$type" '{ print $0,"\t",type }'
}

# in_file=$raw_dir/scap3cd_v1_0.txt.gz
# out_file=$proc_dir/scap3cd_v1_0.txt.gz
# type=$(basename $in_file | cut -d '_' -f 1 | sed 's/scap//')

merged_file=$scap_dir/scap_v1_0.merged.txt.gz

header="chrom\tpos\tref\talt\tscap_score\tscap_type"

# if [[ ! -f $merged_file ]]; then
# 	## Header
# 	echo -e $header | bgzip -c > $merged_file

# 	## Main
# 	counter=0
# 	for in_file in $raw_dir/scap*; do
# 		counter=$((counter+1))
# 		echo "Processing [$counter]: $in_file"
		
# 		type=$(basename $in_file | cut -d '_' -f 1 | sed 's/scap//')

# 		addSpliceType $in_file $type | bgzip -c >> $merged_file
# 	done
# else
# 	echo "SKIPPING: addSpliceType and merge"
# fi

echo -e "\n#========= sort =========#"
job_dir=$base_dir/jobs/; cd $job_dir

sorted_file=$scap_dir/scap_v1_0.sorted.txt.gz

job_file=$job_dir/sortScap.job
echo "module load tabix/1.7" > $job_file
echo "echo -e \"$header\" | bgzip -c > $sorted_file" >> $job_file
echo "zcat $merged_file | tail -n +2 | sort -V -k 1,1 -k 2,2 | bgzip -c >> $sorted_file" >> $job_file

qsub -S /bin/bash -cwd -l h_rt=12:00:00 -l h_vmem=32G -l tmpspace=16G $job_file





