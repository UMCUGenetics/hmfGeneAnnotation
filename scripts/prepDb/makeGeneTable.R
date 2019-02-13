library(biomaRt)
library(readxl)
library(stringr)
options(stringsAsFactors = F)

#========= Paths =========#
base_dir='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'

gene_list <- paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/gene_selection/gene_list.xlsx')
gene_list <- list(
   cancer = as.data.frame(read_excel(gene_list, sheet = 'cancer_genes')),
   hr = as.data.frame(read_excel(gene_list, sheet = 'hr_genes'))
)

hgnc_gene_names <- paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/gene_selection/hgnc_gene_names_20190205.txt')
hgnc_gene_names <- read.delim(hgnc_gene_names,comment.char='#')

mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org')

# ## use this to search for colnames
# findColBiomart <- function(x){
#    att <- listAttributes(mart = mart)
#    att[grep(x,att$name),]
# }
# 
# findColBiomart('hgnc')

#========= Main =========#
getGeneInfo <- function(gene.names){
   #gene.names=gene_list$cancer$gene
   #gene.names=gene_list$hr$gene
   
   hgnc_gene_names <- as.data.frame(apply(hgnc_gene_names,2,toupper))
   gene.names <- toupper(gene.names)
   
   ## Convert gene names to hgnc symbol
   findGene <- function(gene.name, col){
      hgnc_gene_names[which(hgnc_gene_names[col] == gene.name),]
   }
   
   hgnc_gene_info <- do.call(rbind, lapply(gene.names, function(i){
      #i=gene.names[38]
      #i='FANCP'
      
      ## Scan all columns
      df_ss <- findGene(i, 'Approved.symbol')
      if(nrow(df_ss) == 0){ df_ss <- findGene(i, 'Previous.symbol') }
      if(nrow(df_ss) == 0){ df_ss <- findGene(i, 'Alias.symbol') }
      
      ## Output
      if(nrow(df_ss) == 0){
         out <- data.frame(
            gene.name=i,
            hgnc_id=NA,
            hgnc_symbol=NA
         ) 
      } else {
         out <- data.frame(
            gene.name=i,
            hgnc_id=str_extract(df_ss[1,'HGNC.ID'],'\\d+'),
            hgnc_symbol=df_ss[1,'Approved.symbol']
         ) 
      }
      
      return(out)
   }))

   ## Query ensembl
   df <- getBM(
      mart=mart, verbose=F, 
      attributes=c(
         'chromosome_name', 'start_position', 'end_position',
         'hgnc_id','hgnc_symbol','ensembl_gene_id'
      ), 
      filters='hgnc_id', 
      values=hgnc_gene_info$hgnc_id
   )
   
   ## Resolve duplicated. Prioritize smaller entrez gene id
   df <- do.call(rbind, lapply(unique(df$ensembl_gene_id), function(i){
      df_ss <- df[df$ensembl_gene_id == i, ]
      
      if(nrow(df_ss) == 1){ 
         return(df_ss) 
      } else{ 
         return( df_ss[which.min(df_ss$entrezgene),] )
      }
   }))
   
   ## Remove patch entries
   df <- df[grep('PATCH|HSCHR',df$chromosome_name,invert=T),]
   
   ## Reorder
   df <- df[order(df$hgnc_symbol),]
   
   missing_genes <- hgnc_gene_info[!(hgnc_gene_info$hgnc_id %in% df$hgnc_id),'gene.name']
   
   return(list(
      output = df,
      missing.genes = missing_genes
   ))
}

#========= Final gene table =========#
## ignore missing genes
## "RAD51L3-RFFL"  "RP11-286N22.8" "38596"         "40057" 
gene_tables <- list( 
   hr = getGeneInfo(gene_list$hr$gene)$output,
   cancer = getGeneInfo(gene_list$cancer$gene)$output
)

## Remove rows from cancer gene table that are already in hr gene table
gene_tables$cancer <- gene_tables$cancer[
   which( !(gene_tables$cancer$hgnc_symbol %in% gene_tables$hr$hgnc_symbol) )
,]

## Mark as hr gene
gene_tables$hr$is_selected_gene <- TRUE
gene_tables$cancer$is_selected_gene <- FALSE

## Export bed
bed <- do.call(rbind, gene_tables)
colnames(bed)[1] <- paste0('#',colnames(bed)[1])

write.table(
   bed,
   paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/gene_selection/genes.bed'),
   sep='\t',row.names=F,quote=F
)




