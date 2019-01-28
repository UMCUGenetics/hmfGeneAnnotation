library(magrittr)

## Write to clipboard (and paste in excel)
write.clipboard <- function(data){
   clip <- pipe("pbcopy", "w")
   write.table(data, file=clip)
   close(clip)
}

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
   
   # exons <- getBM(mart=mart, verbose=F,
   #                attributes=c("ensembl_gene_id","ensembl_transcript_id","chromosome_name","exon_chrom_start","exon_chrom_end"),
   #                #attributes=c("ensembl_gene_id","ensembl_transcript_id","chromosome_name","3_utr_start","3_utr_end"),
   #                filters='ensembl_gene_id', values=genes$ensembl_gene_id
   # )
   # exon_coords <- merge(genes, exons, by='ensembl_gene_id') %>% .[order(.$hgnc_symbol),]
   # return(exon_coords)
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

#--------- Get all gene coords ---------#
# gene_coords <- getCoordsBiomart()
# gene_coords <- gene_coords[,c('chromosome_name','start_position','end_position','entrezgene','hgnc_symbol','ensembl_gene_id')]
# 
# write.table(gene_coords,
#    '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/data/gene_coords/gene_coords_ensembl.bed',
#    row.names=F, quote=F, sep='\t'
# )

#========= TxDb =========#
# library(GenomicFeatures)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# 
# ## Get HR genes
# getCoordsTxDb <- function(){
#    genome <- TxDb.Hsapiens.UCSC.hg19.knownGene
#    
#    hr_gene_coords <-
#       genes(genome) %>%
#       .[.$gene_id %in% hr_genes$gene_id,] %>%
#       as.data.frame()
#    rownames(hr_gene_coords) <- NULL
#    
#    hr_gene_coords <- merge(hr_genes[,c('gene_id','gene')], hr_gene_coords, by = 'gene_id')
#    
#    ## Get HR gene exons
#    hr_exon_coords <-
#       exons(genome, columns = c('gene_id', 'exon_id')) %>%
#       as.data.frame() %>%
#       .[.$gene_id %in% hr_genes$gene_id,]
#    rownames(hr_exon_coords) <- NULL
#    
#    hr_exon_coords <- merge(hr_genes[,c('gene_id','gene')], hr_exon_coords, by = 'gene_id')
#    return(hr_exon_coords)
# }
# 
# #write.clipboard(hr_gene_coords)
# #write.clipboard(hr_exon_coords)

