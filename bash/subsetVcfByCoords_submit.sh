## submit
hmf_data_dir=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/
manifest_path=${base_dir}/HMF_update/manifest/hmf_file_manifest.txt
variants_dir=${base_dir}/HMF_update/variants_subset/

exist_hr_genes_path=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/metadata/exist_genes/hmf_exist_hr_genes.txt

subsetVcfByCoords_path=${base_dir}/scripts_main/hmfVariantAnnotation/bash/subsetVcfByCoords.sh

job_dir=${base_dir}/HMF_update/jobs/extractVariants/

mkdir -p $job_dir
cd $job_dir

cat $manifest_path | while read sample_name dir multi_vcf_name som_vcf_name gene_cnv_name purity_name; do
	#echo $sample_name

	out_dir=${variants_dir}/$sample_name
	mkdir -p $out_dir

	job_file=${job_dir}/xv_${sample_name}.job
	if [[ -f $job_file ]]; then rm $job_file; fi
	touch $job_file

	## Multi vcf
	multi_vcf_in=${hmf_data_dir}/${dir}/$multi_vcf_name
	multi_vcf_out=${out_dir}/${sample_name}.multi.vcf.gz ## Set multi vcf name to a more regular name

	echo source $subsetVcfByCoords_path > $job_file
	
	echo "echo Processing multi-vcf" >> $job_file
	echo subsetVcfByCoords $multi_vcf_in $multi_vcf_out \"germline\" >> $job_file
	
	## Somatic vcf
	som_vcf_in=${hmf_data_dir}/${dir}/$som_vcf_name
	som_vcf_out=${out_dir}/$som_vcf_name
	
	echo "echo Processing somatic vcf" >> $job_file
	echo subsetVcfByCoords $som_vcf_in $som_vcf_out \"somatic\" >> $job_file

	## Extract cnv for hr genes using list of gene names
	## Gene coords don't completely intersect with those in the purple.gene.cnv file.
	## However, for CNVs precise gene coords are not important.
	cnv_file_in=${hmf_data_dir}/${dir}/$gene_cnv_name
	cnv_file_out=${out_dir}/$gene_cnv_name

	echo "echo Extracting CNVs for $gene_cnv_name" >> $job_file
	## Remove extra colname ExonicBases. 30 colnames, 29 cols + 1 blank col
	echo "cat $cnv_file_in | head -n1 | rev | cut -f2- | rev > $cnv_file_out" >> $job_file
	echo "cat $cnv_file_in | grep -w -f $exist_hr_genes_path >> $cnv_file_out" >> $job_file

	#sh $job_file

	qsub -S /bin/bash -cwd -M N.L.Nguyen-2@umcutrecht.nl -m ea -l h_rt=2:00:00 -l h_vmem=16G $job_file

done