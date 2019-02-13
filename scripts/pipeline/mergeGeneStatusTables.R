args <- commandArgs(trailingOnly=T)

in_dir <- args[1]
out_path <- args[2]
pattern <- ifelse(!is.null(args[3]),args[3],'gene_statuses.txt.gz')


# in_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/'
# pattern <- 'gene_statuses.txt.gz'

files <- system(sprintf("for i in %s/*; do echo $i/%s; done", in_dir, pattern), intern=T)

counter <- 0
n_files <- length(files)

message("\nMerging ",n_files," files...")
pb <- txtProgressBar(max=n_files, style=3)

out <- do.call(rbind, lapply(files, function(i){
   #i=files[1]
   counter <<- counter + 1
   setTxtProgressBar(pb, counter)
   
   sample_name <- gsub('## SAMPLE_NAME: ','',readLines(file(i), n=1))
   sample_name <- strsplit(sample_name,'_')[[1]][2]
   
   df <- read.table(i, sep='\t', header=T)
   
   return(cbind(sample=sample_name, df))
}))

message(sprintf("\nCompressing and writing output..."))
write.table(out, gzfile(out_path), sep='\t', quote=F, row.names=F)
