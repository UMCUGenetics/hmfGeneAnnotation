library(stringr)

## Merge
mergeGeneStatusTables <- function(manifest,out.path,verbose=T){

   #manifest='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/analysis/hr_gene_def/gene_status_paths.txt'
   manifest <- read.table(manifest,sep='\t',header=T,stringsAsFactors=F)
   #manifest <- manifest[1:10,]
   
   #if(dir.exists('/Users/lnguyen/')){ manifest$path <- paste0('/Users/lnguyen/', manifest$path) }
   
   if(verbose){
      counter <- 0
      n_files <- nrow(manifest)
      if(verbose){ message('\nMerging ',n_files,' tables...') }
      pb <- txtProgressBar(max=n_files, style=3)
   }
   
   out <- do.call(rbind, Map(function(path, sample_name){
      counter <<- counter + 1
      #if(verbose){ message('Processing file: ',counter,'/',n_files) }
      #path='/Users/lnguyen//hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD//HMF_update/variants_subset//CPCT02010350R_CPCT02010350T/gene_mut_profile/CPCT02010350T_gene_statuses.txt'
      
      if(verbose){ 
         setTxtProgressBar(pb, counter) 
      }
      
      df <- read.table(path,sep='\t',header=T,stringsAsFactors=F)
      df <- cbind(sample=sample_name, df)
      
   }, manifest$path, manifest$sample_name, USE.NAMES=F))
   
   if(verbose){ message('\nExporting merged table...') }
   #out.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/analysis/hr_gene_def/hmf_gene_statuses.txt'
   write.table(out, out.path, sep='\t', quote=F, row.names=F)
   if(verbose){ message('Done') }
}

args <- commandArgs(trailingOnly=T)
mergeGeneStatusTables(args[1],args[2])
