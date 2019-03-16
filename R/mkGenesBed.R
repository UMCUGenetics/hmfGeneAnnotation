#' Make bed file from a vector of gene names
#'
#' @param gene.names A character vector of gene names
#' @param genes.hgnc A HGNC genes list as a dataframe
#' @param out.path Path to write output. If not provided, a dataframe is returned
#' @param verbose Show messages?
#' @param ... 
#'
#' @return A bed file which also contains the columns: hgnc_id hgnc_symbol ensembl_gene_id
#' @export
#'
#' @example 
#' library(readxl)
#' 
#' base_dir='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'
#' 
#' gene_list <- paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/gene_selection/gene_list.xlsx')
#' gene_list <- list(
#'    cancer = as.data.frame(read_excel(gene_list, sheet = 'cancer_genes')),
#'    selected = as.data.frame(read_excel(gene_list, sheet = 'selected_genes'))
#' )
#' 
#' ## Split selected genes into HR and non-HR genes
#' gene_list$hr <- subset(gene_list$selected, is_hr_gene==T)
#' gene_list$other <- subset(gene_list$selected, is_hr_gene==F)
#' gene_list$selected <- NULL
#' 
#' bed <- do.call(rbind,lapply(gene_list, function(i){
#'    mkGenesBed(i$gene)
#' }))
#' 
#' write.table(
#'    bed,
#'    paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/gene_selection/genes.bed'),
#'    sep='\t',row.names=F,quote=F
#' )


require(biomaRt)
require(stringr)

mkGenesBed <- function(gene.names, genes.hgnc=GENES_HGNC, out.path=NULL, verbose=T, ...){
   #gene.names = gene_list$cancer$gene
   
   if(verbose){ message('Retrieving ENSGs using GENES_HGNC (pre-computed HGNC database) ...') }
   ensg <- geneNamesToEnsg(gene.names, genes.hgnc=genes.hgnc, na.rough.fix = T)
   
   if(verbose){ message('Retrieving gene info using biomaRt...') }
   # ## use this to search for colnames
   # findColBiomart <- function(x){
   #    att <- listAttributes(mart = mart)
   #    att[grep(x,att$name),]
   # }
   # 
   # findColBiomart('hgnc')
   
   mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org',...)
   #mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org')
   
   df <- getBM(
      mart=mart, 
      verbose=F,
      attributes=c(
         'chromosome_name', 'start_position', 'end_position',
         'hgnc_id','hgnc_symbol','ensembl_gene_id'
      ),
      filters='ensembl_gene_id',
      values=ensg
   )
   
   ## Report missing genes
   missing_genes <- ensg[!(ensg %in% df$ensembl_gene_id)]
   if(length(missing_genes) != 0 & verbose){
      warning(sprintf(
         'Info for %s genes could not be retrieved: %s', 
         length(missing_genes), 
         paste(names(missing_genes), collapse=', ')
      ))
   }
   
   ## Remove patch entries
   patch_rows <- grep('PATCH|HSCHR',df$chromosome_name)
   if(length(patch_rows) != 0 & verbose){
      warning(sprintf(
         '%s entries with "PATCH" or "HSCHR" were found. Removing these...', 
         length(patch_rows)
      ))
      df <- df[-patch_rows,]
   }
   
   ## Sort by coordinates
   if(verbose){ message('Sorting genes by coordinates...') }
   chr_order <- unique(str_sort(df$chromosome_name, numeric=T))
   df <- do.call(rbind, lapply(chr_order, function(i){
      df_ss <- df[df$chromosome_name == i,]
      df_ss[order(df_ss$start_position),]
   }))
   
   colnames(df)[1] <- paste0('#',colnames(df)[1])
   
   if(is.null(out.path)){
      return(df)
   } else {
      write.table(df,out.path,sep='\t',row.names=F,quote=F)
   }
}






