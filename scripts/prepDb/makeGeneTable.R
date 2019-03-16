#library(biomaRt)
library(readxl)
library(stringr)
options(stringsAsFactors = F)

#========= Paths =========#
base_dir='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'

gene_list <- paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/gene_selection/gene_list.xlsx')
gene_list <- list(
   cancer = as.data.frame(read_excel(gene_list, sheet = 'cancer_genes')),
   selected = as.data.frame(read_excel(gene_list, sheet = 'selected_genes'))
)

## Split selected genes into HR and non-HR genes
gene_list$hr <- subset(gene_list$selected, is_hr_gene==T)
gene_list$other <- subset(gene_list$selected, is_hr_gene==F)
gene_list$selected <- NULL

hgnc_gene_names <- paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/gene_selection/hgnc_gene_names_20190311.txt')
hgnc_gene_names <- read.delim(hgnc_gene_names,comment.char='#')
colnames(hgnc_gene_names) <- gsub(' ','_',tolower(colnames(GENES_HGNC)))

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
      df_ss <- findGene(i, 'approved_symbol')
      if(nrow(df_ss) == 0){ df_ss <- findGene(i, 'previous_symbols') }
      if(nrow(df_ss) == 0){ df_ss <- findGene(i, 'synonyms') }

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
            hgnc_id=str_extract(df_ss[1,'hgnc_id'],'\\d+'),
            hgnc_symbol=df_ss[1,'approved_symbol']
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
   
   if(length(missing_genes) != 0){
      warning(sprintf(
         '%s missing genes were found: %s', 
         length(missing_genes), 
         paste(missing_genes, collapse=', ')
      ))
   }
   
   return(list(
      output = df,
      missing.genes = missing_genes
   ))
}

#========= Final gene table =========#
## ignore missing genes
## "RAD51L3-RFFL"  "RP11-286N22.8" "38596"         "40057" 

gene_table <- do.call(rbind, lapply(names(gene_list), function(i){
   df <- getGeneInfo(gene_list[[i]]$gene)$output
   df$is_hr_gene <- ifelse(i=='hr',TRUE,FALSE) ## Mark as hr gene
   return(df)
}))

## Remove rows from cancer gene table that are already in gene selection
gene_table <- gene_table[order(gene_table$is_hr_gene, decreasing=T),]
gene_table <- gene_table[!duplicated(gene_table$ensembl_gene_id),]
gene_table <- (function(){
   df <- gene_table
   chr_order <- unique(str_sort(df$chromosome_name, numeric=T))
   do.call(rbind, lapply(chr_order, function(i){
      df_ss <- df[df$chromosome_name == i,]
      df_ss[order(df_ss$start_position),]
   }))
})()

## Export bed
bed <- gene_table
colnames(bed)[1] <- paste0('#',colnames(bed)[1])

write.table(
   bed,
   paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/gene_selection/genes.bed'),
   sep='\t',row.names=F,quote=F
)




