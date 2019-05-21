#!/usr/bin/env Rscript
## Run on hpc
library(devtools)
load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/')

ROOT_DIR <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/'
ini_path_default <- paste0(ROOT_DIR, '/scripts/pipeline/detGeneStatuses_ini.R')

args <- commandArgs(trailingOnly=T)

detGeneStatuses(
   out.dir = args[1],
   cnv.path = args[2],
   germ.path = args[3],
   som.path = args[4],
   purity.path = args[5],
   genes.bed.path = args[6],
   init.path = if(!is.null(args[7])){ args[7] } else { ini_path_default }
)

