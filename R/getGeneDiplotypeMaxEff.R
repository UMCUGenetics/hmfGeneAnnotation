#' Determine the most pathogenic diplotype effect per gene
#'
#' @param gene.diplotypes Table containing the mutation info of alleles 1 and 2
#' @param colname Which column name to take the max of?
#'
#' @return The subsetted gene.diplotypes table with a column indicating how many rows there were
#' originally per gene
#' @export

getGeneDiplotypeMaxEff <- function(gene.diplotypes, colname='hit_score_boosted'){
   #gene.diplotypes=gene_diplotypes
   uniq_genes <- unique(gene.diplotypes$ensembl_gene_id)
   
   l <- lapply(uniq_genes, function(i){
      #print(i)
      #i='ENSG00000139618' ## BRCA2
      df_ss <- gene.diplotypes[gene.diplotypes$ensembl_gene_id==i,]
      hit_scores <- df_ss[,colname]
      
      df_out <- df_ss[which.max(hit_scores),]
      df_out$n_max <- sum(hit_scores==max(hit_scores))
      
      return(df_out)
   })
   
   do.call(rbind,l)
}
