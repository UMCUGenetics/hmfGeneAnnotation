base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/

variants_subset_dir=${base_dir}/HMF_update/variants_subset/
makeVarSigTable_path=${base_dir}/scripts_main/hmfVariantAnnotation/bash/makeVarSigTable.sh

job_dir=${base_dir}/HMF_update/jobs/makeVarSigTables/
mkdir -p $job_dir
cd $job_dir

for i in ${variants_subset_dir}/*/; do
	sample_name=$(basename $i)
	#echo $sample_name

	multi_vcf=${i}/${sample_name}.multi.vcf.gz
	somatic_vcf=${i}/${sample_name}_post_processed_v2.2.vcf.gz

	job_file=${job_dir}/mvst_${sample_name}.job
	touch $job_file

	echo "source $makeVarSigTable_path" > $job_file
	echo "makeVarSigTable $multi_vcf ${i}/germ_varsigs.txt \"germline\"" >> $job_file
	echo "makeVarSigTable $somatic_vcf ${i}/som_varsigs.txt \"somatic\"" >> $job_file

	#qsub -S /bin/bash -cwd -M N.L.Nguyen-2@umcutrecht.nl -m ea -l h_rt=2:00:00 -l h_vmem=16G $job_file
	qsub -S /bin/bash -cwd -m ea -l h_rt=2:00:00 -l h_vmem=16G $job_file
done
