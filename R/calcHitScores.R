#' Sums the pathogenicity score of allele 1 and 2 (hit score)
#'
#' @param gene.diplotypes Table containing the mutation info of alleles 1 and 2
#' @param diplotype.origin.rank A character vector of diplotype origins in descending order
#' of pathogenic impact
#'
#' @return A data frame containing the hit scores and boosted hit scores (scores are boosted by
#' order of diplotype.origin.rank to maintain priority when determining the max diplotype effect)
#' @export
#'
calcHitScores <- function(gene.diplotypes, diplotype.origin.rank=c("cnv_cnv","cnv_som","cnv_germ","germ_som")){
   #gene.diplotypes=gene_diplotypes
   m <- as.matrix(gene.diplotypes[,c('a1.max_score','a2.max_score')])
   hit_score <- unname(rowSums(m))
   
   ## Make boosted hit score for determining the max gene effect (use the 100's place)
   score_boost_lookup <- 100*((length(diplotype.origin.rank)-1):0)
   names(score_boost_lookup) <- diplotype.origin.rank

   score_boost <- score_boost_lookup[
      match(gene.diplotypes$diplotype_origin, names(score_boost_lookup))
   ]

   hit_score_boosted <- hit_score + score_boost
   hit_score_boosted[hit_score==0] <- 0 ## Rm boost for 0 scores

   data.frame(hit_score, hit_score_boosted)
}
