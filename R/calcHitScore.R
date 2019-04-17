#' Calculate hit score that summarizes the biallelic state of a variant pair or gene
#'
#' @param df A dataframe that is the output of mkBialleleMutProfile()
#' @param mode Can be: 'cnv_germ', 'cnv_som', 'germ_som'. Hit score calculation is different for
#' each biallele combination.
#' @param min.hit.score.filter Scores less than this will be set to 0.
#' @param ignore.additional.evidence If TRUE, will only consider the base "5+5" score
#'
#' @return The hit_score as a float
#' @export
#'
calcHitScore <- function(df, mode, min.hit.score.filter=0, ignore.additional.evidence=F){
   #df=biall_mut_profile$cnv_som
   #df=biall_mut_profile$germ_som
   
   modes <- c('cnv_germ','cnv_som','germ_som')
   if(!(mode %in% modes)){
      stop(paste0("Please specify a valid mode: 'cnv_germ', 'cnv_som', 'germ_som'"))
   }
   
   is_def_min_hit_score <- IS_DEF_MIN_HIT_SCORE
   if(ignore.additional.evidence){
      is_def_min_hit_score <- floor(is_def_min_hit_score)
   }
   
   #! NA's have already been set to 0 by mkBialleleMutProfile()
   
   ## Main
   if(mode=='cnv_germ' | mode=='cnv_som'){
      ## Calc hit score
      df$hit_score <- with(df, loh + max_score + alt_exists + ref_loss + cn_break_in_gene)
      df$hit_score[df$full_gene_loss==SCORING$full_gene_loss] <- SCORING$full_gene_loss
      
      ## Determine hit type and def type
      if(mode=='cnv_germ'){ 
         loh_mut_string <- 'loh+germ' 
      } else if(mode=='cnv_som'){ 
         loh_mut_string <- 'loh+som' 
      }
      loh_mut_min_hit_score <- is_def_min_hit_score[[loh_mut_string]]
      
      type_ann <- do.call(rbind, Map(function(full_gene_loss, loh, hit_score){
         def_type <- 'none'
         hit_type <- 'none'
         
         if(full_gene_loss==SCORING$full_gene_loss){
            def_type <- 'full_gene_loss'
            hit_type <- 'full_gene_loss'
         } else {
            if(hit_score >= loh_mut_min_hit_score){ def_type <- loh_mut_string }
            if(loh==SCORING$loh){ hit_type <- loh_mut_string }
         }
         
         return(c(hit_type=hit_type, def_type=def_type))
      }, df$full_gene_loss, df$loh, df$hit_score))
      
      df <- cbind(df, type_ann)
   }
   
   if(mode=='germ_som'){
      df$hit_score <- with(
         df,
         germ.max_score + som.max_score +
         germ.alt_exists + som.alt_exists +
         germ.ref_loss + som.ref_loss +
         cn_break_in_gene
      )
      
      df$hit_type <- 'germ+som'
      
      df$def_type <- 'none'
      df$def_type[df$hit_score >= is_def_min_hit_score['germ+som']] <- 'germ+som' 
   }
   
   ## Set is_def
   df$is_def <- 0
   df$is_def[df$def_type!='none'] <- 1
   
   ## Setting min.hit.score.filter will remove a lot of junk/low impact variants
   df$hit_type[df$hit_score <= min.hit.score.filter] <- 'none'
   
   ## Apply score boosts
   ## Prioritizes higher impact biallelic events when determining most pathogenic biallelic event per gene
   df$hit_score_boosted <- unlist(Map(function(hit_score, def_type){
      if(def_type=='full_gene_loss'){ hit_score + SCORING$score_boost.full_gene_loss }
      else if(def_type=='loh+som'){ hit_score + SCORING$score_boost.loh_som }
      else if(def_type=='loh+germ'){ hit_score + SCORING$score_boost.loh_germ }
      else { hit_score }
   }, df$hit_score, df$def_type))
   
   return(df)
}


# #' Calculate hit score that summarizes the biallelic state of a variant pair or gene
# #'
# #' @param mode Can be: 'cnv_germ', 'cnv_som', 'germ_som'. Hit score calculation is different for
# #' each biallele combination.
# #' @param min.hit.score Scores less than this will be set to 0.
# #' @param full_gene_loss Required by modes containing 'cnv'.
# #' @param loh Required by modes containing 'cnv'.
# #' @param cn_break_in_gene Required by all modes.
# #' @param germ.max_score Required by modes containing 'germ'.
# #' @param germ.alt_exists Required by modes containing 'germ'.
# #' @param germ.ref_loss Required by modes containing 'germ'.
# #' @param som.max_score Required by modes containing 'som'.
# #' @param som.alt_exists Required by modes containing 'som'.
# #'
# #' @return The hit_score as a float
# #' @export
# #'
# calcHitScore <- function(
#    mode=NULL,
#    min.hit.score.filter=0,
#    ignore.additional.evidence=F,
# 
#    ## cnv
#    full_gene_loss=NA, loh=NA, cn_break_in_gene=NA,
# 
#    ## germ
#    germ.max_score=NA, germ.alt_exists=NA, germ.ref_loss=NA,
# 
#    ## som
#    som.max_score=NA, som.alt_exists=NA, som.ref_loss=NA
# ){
#    modes <- c('cnv_germ','cnv_som','germ_som')
#    if(!(mode %in% modes)){
#       stop(paste0("Please specify a valid mode: 'cnv_germ', 'cnv_som', 'germ_som'"))
#    }
# 
#    is_def_min_hit_score <- IS_DEF_MIN_HIT_SCORE
#    if(ignore.additional.evidence){
#       is_def_min_hit_score <- floor(is_def_min_hit_score*10)/10 ## Preserve 1st decimal place (score boosts)
#    }
# 
#    def_type <- 'none'
#    hit_type <- 'none'
# 
#    if(mode=='cnv_germ' | mode=='cnv_som'){
# 
#       if(full_gene_loss==SCORING$full_gene_loss){
#          hit_score <- SCORING$full_gene_loss
#          def_type <- 'full_gene_loss'
#          hit_type <- 'full_gene_loss'
#       }
# 
#       else {
# 
#          if(mode=='cnv_germ'){
#             ## Max clause ensures that germ and som mutations have equal weighting, since
#             ## 'som.ref_loss' does not exist/is not calculated (as it is unreliable), while
#             ## 'germ.ref_loss' does exist
#             # if(anyNA(c(full_gene_loss, loh, germ.max_score, germ.alt_exists, germ.ref_loss, cn_break_in_gene))){
#             #    stop("For mode 'cnv_germ', the required args are: full_gene_loss, loh, germ.max_score, germ.alt_exists, germ.ref_loss, cn_break_in_gene")
#             # }
#             hit_score <- loh + germ.max_score + germ.alt_exists + germ.ref_loss + cn_break_in_gene
#             if(hit_score >= is_def_min_hit_score['loh+germ']){ def_type <- 'loh+germ' }
#             if(loh==SCORING$loh){ hit_type <- 'loh+germ' }
#          }
# 
#          if(mode=='cnv_som'){
#             # if(anyNA(c(full_gene_loss, loh, som.max_score, som.alt_exists, som.ref_loss, cn_break_in_gene))){
#             #    stop("For mode 'cnv_som', the required args are: full_gene_loss, loh, som.max_score, som.alt_exists, cn_break_in_gene")
#             # }
#             hit_score <- loh + som.max_score + som.alt_exists + som.ref_loss + cn_break_in_gene
#             if(hit_score >= is_def_min_hit_score['loh+som']){ def_type <- 'loh+som' }
#             if(loh==SCORING$loh){ hit_type <- 'loh+som' }
#          }
#       }
#    }
# 
#    if(mode=='germ_som'){
#       # if(anyNA(c(germ.max_score, som.max_score,
#       #            germ.alt_exists, som.alt_exists,
#       #            germ.ref_loss, som.ref_loss,
#       #            cn_break_in_gene))
#       # ){
#       #    stop("For mode 'germ_som', the required args are: germ.max_score, som.max_score, germ.alt_exists, som.alt_exists, germ.ref_loss, som.ref_loss, cn_break_in_gene")
#       # }
#       hit_score <-
#          germ.max_score + som.max_score +
#          germ.alt_exists + som.alt_exists +
#          germ.ref_loss + som.ref_loss +
#          cn_break_in_gene
#       if(hit_score >= is_def_min_hit_score['germ+som']){ def_type <- 'germ+som' }
#       hit_type <- 'germ+som'
#    }
# 
#    ## Set is_def
#    is_def <- 0
#    if(def_type!='none'){ is_def <- 1 }
# 
#    ## Setting min.hit.score.filter will remove a lot of junk/low impact variants
#    if(hit_score <= min.hit.score.filter){
#       #hit_score <- 0
#       hit_type <- 'none'
#    }
# 
#    ## Apply score boosts
#    ## Prioritizes higher impact biallelic events when determining most pathogenic biallelic event per gene
#    if(def_type=='full_gene_loss'){ hit_score_boosted <- hit_score + SCORING$score_boost.full_gene_loss }
#    else if(def_type=='loh+som'){ hit_score_boosted <- hit_score + SCORING$score_boost.loh_som }
#    else if(def_type=='loh+germ'){ hit_score_boosted <- hit_score + SCORING$score_boost.loh_germ }
#    else { hit_score_boosted <- hit_score }
# 
#    return(data.frame(
#       hit_score,
#       hit_score_boosted,
#       hit_type,
#       def_type,
#       is_def
#    ))
# }



