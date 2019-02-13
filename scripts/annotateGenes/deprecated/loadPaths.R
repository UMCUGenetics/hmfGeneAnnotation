#========= Set package dir =========#
BASE_DIR <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/'
if(dir.exists('/Users/lnguyen/')){ 
   BASE_DIR <- paste0('/Users/lnguyen/', BASE_DIR) 
}

#========= Dependencies =========#
DEP_DIR <- paste0(BASE_DIR, 'dep/')

DEP_PATHS <- list(
   java = paste0(DEP_DIR,'/jre1.8.0_191/bin/java'),
   snpsift = paste0(DEP_DIR,'/snpEff/SnpSift.jar')
)

#========= Constant paths =========#
FILE_PATHS <- list(
   hr_gene_coords = paste0(BASE_DIR, '/data/gene_coords/hr_gene_coords_ensembl.bed')
)

