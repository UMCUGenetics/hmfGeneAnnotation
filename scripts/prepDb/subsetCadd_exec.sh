#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/
source $ROOT_DIR/loadPaths.sh

job_dir=$ROOT_DIR/jobs/
cd $job_dir

echo "source $(realpath $subsetCadd_sh); subsetCadd $ROOT_DIR/data/variant_significance/cadd/InDels_inclAnno.tsv.gz" > $job_dir/ss_cad_indel.job
echo "source $(realpath $subsetCadd_sh); subsetCadd $ROOT_DIR/data/variant_significance/cadd/whole_genome_SNVs_inclAnno.tsv.gz" > $job_dir/ss_cad_snv.job

#qsub -S /bin/bash -cwd -l h_rt=24:00:00 -l h_vmem=16G $job_dir/ss_cad_indel.job
qsub -S /bin/bash -cwd -l h_rt=24:00:00 -l h_vmem=50G -pe threaded 24 $job_dir/ss_cad_snv.job