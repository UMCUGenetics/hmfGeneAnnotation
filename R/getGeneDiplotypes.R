#' Per gene, get row(s) where score is max
#'
#' @param df A dataframe containing ensembl_gene_id and the desired column for subsetting
#' @param mode Can be: 'cnv_germ', 'cnv_som', 'germ_som'.
#' @param simplify.snpeff.eff Simplify snpeff_eff (e.g. stop_gained, stop_loss; becomes nonsense)
#'
#' @return A subsetting dataframe
#' @export
#'
getGeneDiplotypes <- function(df, mode, simplify.snpeff.eff=T){
   
   #df=biall_mut_profile$cnv_som
   
   SEL_COLS <- list(
      common=c('ensembl_gene_id','hgnc_symbol','hit_score','def_type','is_def'),
      allele=c('chrom','pos','hgvs_c','max_score','max_score_origin')
   )

   modes <- c('cnv_germ','cnv_som','germ_som')
   if(!(mode %in% modes)){
      stop(paste0("Please specify a valid mode: 'cnv_germ', 'cnv_som', 'germ_som'"))
   }

   simplifySnpeffEff <- function(v){
      v_new <- SNPEFF_SIMPLE_ANN_LOOKUP[match(v, SNPEFF_SIMPLE_ANN_LOOKUP$ann), 'ann_s2']
      v_new[is.na(v_new)] <- 0
      return(v_new)
   }

   if(mode=='cnv_germ' | mode=='cnv_som'){
      diplotypes <- do.call(rbind, lapply(1:nrow(df), function(i){
         row <- df[i,]
         if(row$full_gene_loss==SCORING$full_gene_loss){
            a1 <- a2 <- 'full_gene_loss'
            a1.origin <- a2.origin <- 'cnv'
         } else if(row$loh==SCORING$loh) {
            a1 <- 'loh'
            a1.origin <- 'cnv'
            a2 <- row$snpeff_eff
            a2.origin <- if(mode=='cnv_germ'){ 'germ' } else if(mode=='cnv_som'){ 'som' }
         } else {
            a1 <- 'none'
            a1.origin <- 'cnv'
            a2 <- row$snpeff_eff
            a2.origin <- if(mode=='cnv_germ'){ 'germ' } else if(mode=='cnv_som'){ 'som' }
         }

         return(data.frame(a1,a2,a1.origin,a2.origin))
      }))


      out_a1 <- (function(){
         col_names <- c(paste0('a1.', SEL_COLS$allele))

         out_a1 <- data.frame(matrix(nrow=nrow(df),ncol=length(col_names)))
         colnames(out_a1) <- col_names

         out_a1 <- cbind(out_a1, diplotypes[,c('a1','a1.origin')])

         return(out_a1)
      })()

      out_a2 <- df[,SEL_COLS$allele]
      colnames(out_a2) <- paste0('a2.',colnames(out_a2))
      
      out_a2 <- cbind(out_a2, diplotypes[,c('a2','a2.origin')])
      
      out_a2$a2 <- 
         if(simplify.snpeff.eff){
            simplifySnpeffEff(diplotypes$a2)
         } else {
            diplotypes$a2
         }

      out <- cbind(df[,SEL_COLS$common], diplotype_origin = mode, out_a1, out_a2)
   }

   if(mode=='germ_som'){
      out <- cbind(
         df[,SEL_COLS$common],
         diplotype_origin = mode,
         
         df[,paste0('germ.',SEL_COLS$allele)],
         a1=df$germ.snpeff_eff,
         a1.origin='germ',
         
         df[,paste0('som.',SEL_COLS$allele)],
         a2=df$som.snpeff_eff,
         a2.origin='som'
      )

      colnames(out) <- gsub('germ','a1',colnames(out))
      colnames(out) <- gsub('som','a2',colnames(out))

      if(simplify.snpeff.eff){
         out$a1 <- simplifySnpeffEff(out$a1)
         out$a2 <- simplifySnpeffEff(out$a2)
      }
   }

   return(out)
}

# getVarPairDiplotype <- function(
#    mode,
#    full_gene_loss=NA, loh=NA,
#    a1.snpeff_eff=NA, a1.max_score=NA,
#    a2.snpeff_eff=NA, a2.max_score=NA
# ){
#    modes <- c('cnv_germ','cnv_som','germ_som')
#    if(!(mode %in% modes)){
#       stop(paste0("Please specify a valid mode: 'cnv_germ', 'cnv_som', 'germ_som'"))
#    }
#
#    if(mode=='cnv_germ' | mode=='cnv_som'){
#       if(anyNA(c(full_gene_loss, loh, a2.max_score, a2.snpeff_eff))){
#          stop("For modes 'cnv_germ' or 'cnv_som', the required args are: full_gene_loss, loh, a2.max_score, a2.snpeff_eff")
#       }
#       if(full_gene_loss==SCORING$full_gene_loss){
#          a1 <- a2 <- 'full_gene_loss'
#       } else if(loh==SCORING$loh) {
#          a1 <- 'loh'
#          a2 <- a2.snpeff_eff
#       } else {
#          a1 <- 'none'
#          a2 <- a2.snpeff_eff
#       }
#    }
#
#    if(mode=='germ_som'){
#       if(anyNA(c(a1.snpeff_eff, a2.snpeff_eff, a1.max_score, a2.max_score))){
#          stop("For modes 'cnv_germ' or 'cnv_som', the required args are: a1.max_score, a2.max_score, a1.snpeff_eff, a2.snpeff_eff")
#       }
#       # max_scores <- c(a1.max_score, a2.max_score)
#       # snp_effs <- c(a1.snpeff_eff, a2.snpeff_eff)
#       #
#       # a1 <- snp_effs[which.max(max_scores)]
#       # a2 <- snp_effs[which.min(max_scores)]
#
#       a1 <- a1.snpeff_eff
#       a2 <- a2.snpeff_eff
#    }
#
#    return(data.frame(
#       a1, a2, a1.max_score, a2.max_score
#    ))
# }
