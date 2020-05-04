####################################################################################################
#' Add annotations to the gene cnv table
#'
#' @description Annotations deep deletions, truncations, LOH, and amplifications
#'
#' @param gene.cnv purple gene cnv dataframe
#' @param arm.ploidies arm ploidies dataframe
#' @param min.cn.arm.ploidy.diff Minimum difference between min_copy_number and arm_ploidy before a
#' gain is considered to have happened
#' @param deep.del.max.max.copy.number The max max_copy_number for a gene to be considered to 
#' be completely lost (deep deletion)
#' @param trunc.max.min.copy.number The max min_copy_number to for a gene to be considered to have
#' a truncation
#' @param loh.max.min.minor.allele.ploidy The maximum min_minor_allele_ploidy for a gene to be
#' considered to have loss of heterozygosity
#' @param scoring a list in the form list(deep_deletion=c(5,5),trunc=c(5,5),loh=c(5,0))
#' @param trunc.is.deep.del Consider truncations the same as deep deletions?
#'
#' @return The original input dataframe with annotation columns
#' @export
#'
mkMutProfileGeneCnv <- function(
   gene.cnv, arm.ploidies,
   
   min.cn.arm.ploidy.diff=0.8,
   
   ## Cutoffs
   deep.del.max.max.copy.number = 0.3, 
   trunc.max.min.copy.number = 0.3, 
   loh.max.min.minor.allele.ploidy = 0.2,

   ## Scoring
   scoring=list(deep_deletion=c(5,5),trunc=c(5,5),loh=c(5,0)),
   trunc.is.deep.del=T,
   
   ## Misc
   verbose=T
){
   if(nrow(gene.cnv)==0){ stop('No rows are present in gene cnv table') }
   if(nrow(arm.ploidies)==0){ stop('No rows are present in arm ploidies table') }

   #gene.cnv=input_tables$gene_cnv
   #arm.ploidies=input_tables$arm_ploidies
   
   #--------- Gains ---------#
   if(verbose){ message('> Gains...') }
   gains <- data.frame(
      min_copy_number = gene.cnv$min_copy_number,
      genome_ploidy = arm.ploidies[arm.ploidies$chrom=='genome','ploidy'],
      arm_ploidy = arm.ploidies$ploidy[ match(gene.cnv$chrom_arm, arm.ploidies$chrom) ],
      stringsAsFactors = F
   )
   
   if(verbose){ message('Determining amplification type...') }
   gains$cn_arm_ploidy_diff <- abs(gains$min_copy_number - gains$arm_ploidy)
   
   gains$amp_type <- 'none'
   gains <- within(gains,{
      amp_type[ arm_ploidy > genome_ploidy ] <- 'arm'
      amp_type[ min_copy_number > arm_ploidy & cn_arm_ploidy_diff >= min.cn.arm.ploidy.diff ] <- 'local'
      #amp_type[ amp_type=='local' & min_copy_number < genome_ploidy ] <- 'local.lt.genome'
   })
   
   if(verbose){ message('Calculating amplification level...') }
   gains <- within(gains,{
      amp_ratio_arm <- arm_ploidy / genome_ploidy
      amp_ratio_arm[amp_type!='arm'] <- 1
      
      amp_ratio_local <- min_copy_number / arm_ploidy
      amp_ratio_local[amp_type!='local'] <- 1
   })
   
   gains$amp_ratio <- unlist(Map(function(arm, local, amp_type){
      if(amp_type=='arm'){ return(arm) }
      if(amp_type=='local'){ return(local) }
      return(1)
   }, gains$amp_ratio_arm, gains$amp_ratio_local, gains$amp_type, USE.NAMES=F))
   
   # if(verbose){ message('Assigning amplification bin...') }
   # gain_level_names <- paste0(
   #    '_[',amp.bins[-length(amp.bins)],
   #    ',',amp.bins[-1],')'
   # )
   # gain_level_names[ .bincode(2, amp.bins, include.lowest=T, right=F) ]
   # 
   # gains$amp_context <- gain_level_names[ .bincode(gains$amp_ratio, amp.bins, include.lowest=T) ]
   # gains$amp_context[gains$amp_type=='none'] <- ''
   # gains$amp_context <- paste0(gains$amp_type,gains$amp_context)
   
   #--------- Losses ---------#
   if(verbose){ message('> Losses...') }
   losses <- data.frame(matrix(0,nrow=nrow(gene.cnv), ncol=3), stringsAsFactors=F)
   colnames(losses) <- c('deep_deletion','trunc','loh')
   
   if(verbose){ message('Determining the presence of inactivating CNV events...') }
   losses <- within(losses,{
      deep_deletion[ gene.cnv$max_copy_number <= deep.del.max.max.copy.number ] <- 1
      trunc[ gene.cnv$min_copy_number <= trunc.max.min.copy.number ] <- 1
      loh[ gene.cnv$min_minor_allele_ploidy <= loh.max.min.minor.allele.ploidy ] <- 1
   })
   
   if(verbose){ message('Determining most impactful loss type...') }
   losses$loss_type <- (function(){
      v <- colnames(losses)[ max.col(losses, ties.method='first') ]
      v[rowSums(losses)==0] <- 'none'
      return(v)
   })()
   
   if(trunc.is.deep.del){
      losses$loss_type[losses$loss_type=='trunc'] <- 'deep_deletion'
   }
   
   if(verbose){ message('Assigning biallele score to loss events...') }
   loss_scores <- unname( scoring[match(losses$loss_type, names(scoring))] )
   loss_scores <- do.call(rbind, lapply(loss_scores, function(i){
      if(is.null(i)){ c(0,0) }
      else { i }
   }))
   colnames(loss_scores) <- c('a1.score','a2.score')
   
   out <- cbind(
      gene.cnv, 
      subset(gains,select=-min_copy_number),
      cn_loss_type=losses$loss_type,loss_scores
   )
   
   ## Shorten floating point numbers
   out <- as.data.frame(lapply(out, function(i){
      if(is.double(i)){ round(i,3) }
      else{ i }
   }))
   
   return(out)
}

####################################################################################################
#' Create annotations to the germline/somatic tables originating from HMF germline/somatic vcfs
#'
#' @description Primarily used to determine
#'
#' @param df.snv.indel A germline/somatic table originating from HMF germline/somatic vcfs
#' @param scoring A list containing the scoring for snpeff, clinvar, and enigma annotation values
#' @param keep.only.first.eff Only keep first snpeff_eff (items are separated by '&')?
#' @param verbose Show messages?
#'
#' @return The original input dataframe with annotation columns
#' @export
#'
mkMutProfileSnvIndel <- function(
   df.snv.indel, 
   scoring=SCORING_MUT,
   keep.only.first.eff=T,
   filter.no.impact.variants=T,
   verbose=T
){
   
   #--------- Sanity checks ---------#
   if(nrow(df.snv.indel)==0){ 
      return(data.frame()) 
   }
   #df.snv.indel=input_tables$som_txt
   #df.snv.indel=input_tables$germ_txt
   
   #--------- Score annotations ---------#
   if(verbose){ message('Getting snpEff scores...') }
   snpeff_eff <- sapply(strsplit(df.snv.indel$snpeff_eff,'&'),`[[`,1)
   
   snpeff_score <- unname(scoring$snpeff)[ match(snpeff_eff,names(scoring$snpeff)) ]
   snpeff_score[is.na(snpeff_score)] <- 0 ## If annotation not found in snpeff table, return 0
   
   if(verbose){ message('Getting ClinVar scores...') }
   clinvar_score <- unname(scoring$clinvar)[ match(df.snv.indel$clinvar_sig,names(scoring$clinvar)) ]
   clinvar_score[is.na(clinvar_score)] <- 0
   
   #--------- Calculate max score and which database it came from ---------#
   if(verbose){ message('Calculating max scores...') }
   sig_scores <- data.frame(clinvar_score, snpeff_score, stringsAsFactors=F)
   
   sig_scores$max_score <- unlist(Map(function(clinvar_score, snpeff_score){
      if(clinvar_score != 0){ return(clinvar_score) }
      if(snpeff_score != 0){ return(snpeff_score) }
      return(0)
   }, clinvar_score, snpeff_score, USE.NAMES=F))
   
   sig_scores$max_score_origin <- unlist(Map(function(clinvar_score, snpeff_score){
      if(clinvar_score != 0){ return('clinvar') }
      if(snpeff_score != 0){ return('snpeff') }
      return('none')
   }, clinvar_score, snpeff_score, USE.NAMES=F))
   
   #--------- Output ---------#
   out <- cbind(df.snv.indel, sig_scores)
   
   ## Deal with multiple sneff effs
   if(keep.only.first.eff){
      if(verbose){ message('Keeping only first snpeff_eff...') }
      out$snpeff_eff <- gsub('&.+$','',out$snpeff_eff)
   }
   
   if(filter.no.impact.variants){
      out$snpeff_eff[ 
         !(out$snpeff_eff %in% names(scoring$snpeff)) & out$max_score==0
      ] <- 'no_impact_variant'
   }
   
   return(out)
}



