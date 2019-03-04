options(stringsAsFactors=F)

args <- commandArgs(trailingOnly=TRUE)

gene_cnv_path <- args[1]
out_path <- args[2]

ROOT_DIR <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/'
GENES_ENST2ENSG <- paste0(ROOT_DIR,'/data/gene_selection/human_genes_enst2ensg.txt.gz')
GENES_BED <- paste0(ROOT_DIR,'/data/gene_selection/genes.bed')

# gene_cnv_path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/170416_HMFregCPCT_FR12245150_FR14064764_CPCT02020493/CPCT02020493T.purple.gene.cnv'
# GENES_ENST2ENSG='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/data/gene_selection/human_genes_enst2ensg.txt.gz'
# GENES_BED='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/data/gene_selection/genes.bed'

gene_cnv <- read.delim(gene_cnv_path)
human_genes_enst2ensg <- read.delim(gzfile(GENES_ENST2ENSG), check.names=F)
genes_bed <- read.delim(GENES_BED, check.names=F)

## Convert colnames from camelcase to snakecase
colnames(gene_cnv) <- gsub("([a-z])([A-Z])", "\\1_\\L\\2", colnames(gene_cnv), perl = TRUE)
colnames(gene_cnv) <- tolower(colnames(gene_cnv))

## Retrieve ensembl gene ids
message('Retrieving ensembl gene ids...')
# gene_info_ensembl <- getBM(
#    mart=useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org'),
#    attributes=c('hgnc_symbol','ensembl_gene_id','ensembl_transcript_id'),
#    filters='ensembl_transcript_id', 
#    values=gene_cnv$transcript_id,
#    verbose=F
# )

ins_cols <- human_genes_enst2ensg[
   match(gene_cnv$transcript_id, human_genes_enst2ensg$ensembl_transcript_id),
   c('ensembl_gene_id','hgnc_symbol')
]
gene_cnv <- cbind(gene_cnv,ins_cols)

## Subset by ensembl gene id of selected genes
message('Subsetting gene cnv table...')
gene_cnv <- gene_cnv[gene_cnv$ensembl_gene_id %in% genes_bed$ensembl_gene_id,]

## Remove duplicate gene ids
gene_cnv <- gene_cnv[!(duplicated(gene_cnv$ensembl_gene_id)),]

## Select relevant rows
sel_cols <- c(
   'chromosome','start','end','hgnc_symbol','ensembl_gene_id',
   'min_copy_number','max_copy_number','min_minor_allele_ploidy'
)
gene_cnv <- gene_cnv[,sel_cols]
colnames(gene_cnv)[1] <- 'chrom'

## Export
message('Exporting new gene cnv table...')
write.table(gene_cnv, out_path, sep='\t', row.names=F, quote=F)


