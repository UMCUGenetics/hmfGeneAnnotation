#' Calculate hit score that summarizes the biallelic state of a variant pair or gene
#'
#' @param mode Can be: 'cnv_germ', 'cnv_som', 'germ_som'. Hit score calculation is different for
#' each biallele combination.
#' @param min.hit.score Scores less than this will be set to 0.
#' @param full_gene_loss Required by modes containing 'cnv'.
#' @param loh Required by modes containing 'cnv'.
#' @param cn_break_in_gene Required by all modes.
#' @param germ.max_score Required by modes containing 'germ'.
#' @param germ.alt_exists Required by modes containing 'germ'.
#' @param germ.ref_loss Required by modes containing 'germ'.
#' @param som.max_score Required by modes containing 'som'.
#' @param som.alt_exists Required by modes containing 'som'.
#'
#' @return The hit_score as a float
#' @export
#'

calcHitScore <- function(
   mode=NULL,
   min.hit.score=0,
   ignore.additional.evidence=F,

   ## cnv
   full_gene_loss=NA, loh=NA, cn_break_in_gene=NA,

   ## germ
   germ.max_score=NA, germ.alt_exists=NA, germ.ref_loss=NA,

   ## som
   som.max_score=NA, som.alt_exists=NA
){
   modes <- c('cnv_germ','cnv_som','germ_som')
   if(!(mode %in% modes)){
      stop(paste0("Please specify a valid mode: 'cnv_germ', 'cnv_som', 'germ_som'"))
   }

   is_def_min_hit_score <- IS_DEF_MIN_HIT_SCORE
   if(ignore.additional.evidence){
      is_def_min_hit_score <- floor(is_def_min_hit_score*10)/10 ## Preserve 1st decimal place (score boosts)
   }

   def_type <- 'none'
   if(mode=='cnv_germ' | mode=='cnv_som'){

      if(full_gene_loss==SCORING$full_gene_loss){
         hit_score <- SCORING$full_gene_loss
         def_type <- 'full_gene_loss'
      }

      else {

         if(mode=='cnv_germ'){
            ## Max clause ensures that germ and som mutations have equal weighting, since
            ## 'som.ref_loss' does not exist/is not calculated (as it is unreliable), while
            ## 'germ.ref_loss' does exist
            if(anyNA(c(full_gene_loss, loh, germ.max_score, germ.alt_exists, cn_break_in_gene, germ.ref_loss))){
               stop("For mode 'cnv_germ', the required args are: full_gene_loss, loh, germ.max_score, germ.alt_exists, germ.ref_loss, cn_break_in_gene")
            }
            hit_score <- loh + germ.max_score + germ.alt_exists + max(germ.ref_loss, cn_break_in_gene)
            if(hit_score >= is_def_min_hit_score['loh+germ']){ def_type <- 'loh+germ' }
         }

         if(mode=='cnv_som'){
            if(anyNA(c(full_gene_loss, loh, som.max_score, som.alt_exists, cn_break_in_gene))){
               stop("For mode 'cnv_som', the required args are: full_gene_loss, loh, som.max_score, som.alt_exists, cn_break_in_gene")
            }
            hit_score <- loh + som.max_score + som.alt_exists + cn_break_in_gene
            if(hit_score >= is_def_min_hit_score['loh+som']){ def_type <- 'loh+som' }
         }
      }
   }

   if(mode=='germ_som'){
      if(anyNA(c(germ.max_score, som.max_score, germ.alt_exists, som.alt_exists, cn_break_in_gene, germ.ref_loss))){
         stop("For mode 'germ_som', the required args are: germ.max_score, som.max_score, germ.alt_exists, som.alt_exists, cn_break_in_gene, germ.ref_loss")
      }
      hit_score <-
         germ.max_score + som.max_score +
         germ.alt_exists + som.alt_exists +
         cn_break_in_gene + germ.ref_loss
      if(hit_score >= is_def_min_hit_score['germ+som']){ def_type <- 'germ+som' }
   }

   ## Apply score boosts
   if(def_type=='full_gene_loss'){ hit_score <- hit_score + SCORING$score_boost.full_gene_loss }
   if(def_type=='loh+som'){ hit_score <- hit_score + SCORING$score_boost.loh_som }
   if(def_type=='loh+germ'){ hit_score <- hit_score + SCORING$score_boost.loh_germ }

   ## Other post processing
   is_def <- 0
   if(def_type!='none'){ is_def <- 1 }

   if(hit_score==SCORING$loh){ def_type <- 'loh_only' }
   if(hit_score < min.hit.score){ hit_score <- 0 }

   return(data.frame(
      hit_score,
      def_type,
      is_def
   ))
}



