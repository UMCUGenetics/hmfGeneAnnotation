base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/

manifest_path=${base_dir}/HMF_update/manifest/hmf_file_manifest.txt
variants_subset_dir=${base_dir}/HMF_update/variants_subset/
hmf_data_dir=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/

cat $manifest_path | while read sample_name dir_name multi_vcf som_vcf gene_cnv purity; do
	echo $sample_name
	cp ${hmf_data_dir}/${dir_name}/$purity ${variants_subset_dir}/${sample_name}/$purity
done
