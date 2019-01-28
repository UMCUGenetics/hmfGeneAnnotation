library(stringr)

purity_files <- list.files('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/', 
                           pattern = '.purple.purity$', recursive = T, full.names = T)

purity <- unlist(lapply(purity_files, function(i){ 
   message(i)
   read.table(i, header = T, comment.char='', check.names = F)[[1]]
}))

out <- data.frame(
   sample_id = str_remove(basename(purity_files), '.purple.purity'),
   tumor_purity = purity
)

write.table(out, '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/metadata/hrd_gene_annotation/tumor_purity.txt')
