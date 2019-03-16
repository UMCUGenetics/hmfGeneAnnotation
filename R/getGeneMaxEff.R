#' Return row(s) with the maximum effect score per gene
#'
#' @param df A mut_profile dataframe containing a max_score column
#' @param colname An colname to determine the max effect. Default='max_score'
#' @param greedy.max.eff Return only the one instance with maximum score per gene?
#' @param show.n.max Add a column indicating the number of rows with max score per gene?
#'
#' @return A subsetted dataframe
#' @export
#'
getGeneMaxEff <- function(df, colname='max_score', greedy.max.eff=T, show.n.max=F){
   #df=gene_mut_profile$germ
   if(nrow(df) == 0){
      stop(return(NA))
   }

   gene_ids <- unique(df$ensembl_gene_id)
   do.call(rbind, lapply(gene_ids, function(i){
      #i="BRE"

      df_ss <- df[ df$ensembl_gene_id == i,]

      max_max_score <- max(df_ss[,colname])
      df_ss <- df_ss[df_ss[,colname] == max_max_score, ]


      if(greedy.max.eff){
         out <- df_ss[1,]

         if(show.n.max){
            n_max <- nrow(df_ss)
            out$n_max <- n_max
         }
      }

      return(out)
   }))
}
