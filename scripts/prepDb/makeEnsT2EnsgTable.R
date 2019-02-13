library(biomaRt)
options(stringsAsFactors = F)

base_dir='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'

mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org')
enst2ensg <- getBM(
   mart=mart, verbose=F, 
   attributes=c('hgnc_symbol','ensembl_gene_id','ensembl_transcript_id')
)

write.table(
   enst2ensg,
   gzfile(paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/gene_selection/human_genes_enst2ensg.txt.gz')),
   sep='\t',row.names=F,quote=F
)