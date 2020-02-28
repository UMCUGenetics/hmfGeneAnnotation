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
calcHitScores <- function(
   gene.diplotypes, 
   diplotype.origin.rank=c('cnv_cnv','cnv_som','cnv_germ','som_som','germ_som')
){
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

####################################################################################################
#' Determine biallelic hit type
#'
#' @param diplotypes Diplotypes dataframe
#' @param min.a1.score Min score to assess hit type
#' @param min.a2.score See above
#'
#' @return A character vector of hit types
#' @export
#'
detHitType <- function(diplotypes, min.a1.score=1, min.a2.score=1){
   
   diplotypes_ss <- as.matrix(
      diplotypes[ c('a1.eff','a2.eff','a1.max_score','a2.max_score','diplotype_origin') ]
   )
   rownames(diplotypes_ss) <- NULL
   
   apply(diplotypes_ss, 1, function(i){
      a1.eff <- i[1]
      a2.eff <- i[2]
      a1.max_score <- i[3]
      a2.max_score <- i[4]
      diplotype_origin <- i[5]
      
      if(a1.eff %in% c('deep_deletion','trunc')){
         return(a1.eff)
      } 
      
      if(a1.eff=='loh'){
         if(diplotype_origin=='cnv_germ' & a2.max_score>=min.a2.score){ return('loh+germ') }
         if(diplotype_origin=='cnv_som' & a2.max_score>=min.a2.score){ return('loh+som') }
         return('loh_only')
      } 
      
      if(
         diplotype_origin %in% c('germ_som','som_som') &
         a1.max_score>=min.a1.score & a2.max_score >=min.a2.score
      ){
         if(diplotype_origin=='germ_som'){ return('germ+som') }
         if(diplotype_origin=='som_som'){ return('som+som') }
      } 
      
      return('none')
   })
}

####################################################################################################
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


####################################################################################################
#' Determines the reading frame in case of one or more frameshifts
#'
#' @param gene.diplotypes Table containing the mutation info of alleles 1 and 2
#' @param mut.profile.germ Germ mut profile table
#' @param mut.profile.som Som mut profile table
#'
#' @return A dataframe containing info of: the number of detected indels, the lengths of each
#' frameshift, the reading frame (sum of the frameshift lengths), and a logical indicating whether
#' the gene is out of frame
#' @export
#'
detReadingFrame <- function(gene.diplotypes, mut.profile.germ, mut.profile.som){
   # gene.diplotypes=gene_diplotypes_max
   # mut.profile.germ=mut_profile$germ
   # mut.profile.som=mut_profile$som
   
   has_fs <- with(gene.diplotypes,{ a1.eff=='frameshift_variant' | a2.eff=='frameshift_variant'})
   fs_genes <- gene.diplotypes[has_fs,'hgnc_symbol']
   
   ## Only look back into mut_profile if gene has frameshift
   fs_info <- do.call(rbind,lapply(fs_genes, function(i){
      #i='BRCA2'
      mut_profile_ss <- lapply(list(mut.profile.germ, mut.profile.som), function(j){
         j[with(j,{ hgnc_symbol==i & snpeff_eff=='frameshift_variant'}),]
      })
      # mut_profile_ss[[1]]$mut_origin <- 'germ'
      # mut_profile_ss[[2]]$mut_origin <- 'som'
      
      df <- rbind(mut_profile_ss[[1]], mut_profile_ss[[2]])
      
      ##
      fs_lengths <- nchar(df$alt) - nchar(df$ref)
      reading_frame <- sum(fs_lengths)
      if(reading_frame %% 3 == 0){
         out_of_frame <- 0
      } else {
         out_of_frame <- 1
      }
      
      ##
      fs_origins <- c(
         rep('germ',nrow(mut_profile_ss[[1]])),
         rep('som',nrow(mut_profile_ss[[2]]))
      )
      
      ## 
      fs_positions <- c(
         mut_profile_ss[[1]]$pos,
         mut_profile_ss[[2]]$pos
      )
      
      return(data.frame(
         n_fs=nrow(df),
         fs_lengths=paste(fs_lengths,collapse=','),
         fs_origins=paste(fs_origins, collapse=','),
         fs_positions=paste(fs_positions,collapse=','),
         reading_frame,
         out_of_frame
      ))
   }))
   
   ## Initialize output
   out <- data.frame(
      n_frameshifts=0,
      fs_lengths='',
      fs_origins='',
      fs_positions='',
      reading_frame=0,
      out_of_frame=0
   )
   out <- out[rep(1,nrow(gene.diplotypes)),]
   
   ## Fill in rows (genes) which have frameshift
   out[has_fs,] <- fs_info
   
   return(out)
}
