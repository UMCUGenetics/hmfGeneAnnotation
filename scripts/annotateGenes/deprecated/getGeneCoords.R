library(magrittr)

## Load list of hr genes
hr_genes <- read.table(
   '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/data/hr_genes.txt',
   sep='\t', header = T, check.names = F, stringsAsFactors = F
)

#========= ENSEMBL =========#
## Use ENSEMBL coordinates. SnfEff in HMF pipeline uses ENSEMBL data for annotations
library(biomaRt)

mart <- useMart(
   biomart = "ensembl", dataset = "hsapiens_gene_ensembl", 
   host="grch37.ensembl.org"
)

## use this to search for colnames
findColBiomart <- function(x){
   att <- listAttributes(mart = mart)
   att[grep(x,att$name),]
}

#findColBiomart('start_position')

getCoordsBiomart <- function(entrez.id = NULL){
   
   params <- list(
      attributes = c("entrezgene","hgnc_symbol","ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
      filters = 'entrezgene',
      values=entrez.id
   )
   
   if(is.null(entrez.id)){
      genes <- getBM(mart=mart, verbose=F,attributes=params$attributes)
   } else {
      genes <- getBM(
         mart=mart, verbose=F, attributes=params$attributes, 
         filters=params$filters, values=params$values
      )
   }
   
   genes <- genes %>% 
      #.[!is.na(.$entrezgene),] %>%
      .[.$chromosome_name %in% c(1:22, 'X','Y'),] %>% ## keep only autosomal and sex genes
      .[nchar(.$hgnc_symbol) != 0,] %>% ## keep only genes with gene symbols (likely won't investigate others)
      .[order(.$hgnc_symbol),] %>% ## sort by gene name
      .[!grepl('PATCH',.$chromosome_name),] ## remove rows where CHROM is 'PATCH'
   
   return(genes)
}

#--------- Get hr gene coords ---------#
## Using ENSEMBL exon coords still gives intron_variant in exons
## Better to just include all variants in a gene region
gene_coords <- getCoordsBiomart(hr_genes$gene_id)
gene_coords <- gene_coords[,c('chromosome_name','start_position','end_position','entrezgene','hgnc_symbol','ensembl_gene_id')]
colnames(gene_coords)[1] <- '#chromosome_name'

## Add aliases
gene_coords$aliases <- unlist(lapply(gene_coords$entrezgene, function(i){
   hr_genes[hr_genes$gene_id == i,'string']
}))

## Export
write.table(gene_coords,
   '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/data/gene_coords/hr_gene_coords_ensembl.bed',
   row.names=F, quote=F, sep='\t'
)