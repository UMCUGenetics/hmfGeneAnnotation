#' Create annotations to the germline/somatic tables originating from HMF germline/somatic vcfs
#'
#' @description Primarily used to determine
#'
#' @param df.snv.indel A germline/somatic table originating from HMF germline/somatic vcfs
#' @param tumor.purity Tumor purity
#' @param mode Can be 'germline' or 'somatic'
#' @param rm.non.selected.genes Remove genes not specified by the user in genes.bed
#' @param genes.bed A bed file which contains ENSEMBL gene ids
#' @param gene.identifier Can be 'snpeff_gene' or 'ensembl_gene_id'. If 'snpeff_gene', then ENSEMBL
#' gene ids will be retrieved
#' @param keep.only.first.eff Only keep first snpeff_eff (items are separated by '&')?
#' @param verbose Show messages?
#'
#' @return The original input dataframe with annotation columns
#' @export
#'

# sample_name='CPCT02070023R_CPCT02070023TII'
# out_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/',sample_name)
#
# input_paths <- list(
#    cnv = paste0(out_dir,'/',sample_name,'.purple.gene.cnv'),
#    germ = paste0(out_dir,'/varsig/',sample_name,'_varsigs_germ.txt.gz'),
#    som = paste0(out_dir,'/varsig/',sample_name,'_varsigs_som.txt.gz'),
#    purity = paste0(out_dir,'/',sample_name,'.purple.purity')
# )
#
# input <- list(
#    cnv = read.delim(input_paths$cnv),
#    germ = read.delim(input_paths$germ),
#    som = read.delim(input_paths$som),
#    purity = read.table(input_paths$purity, skip=1)[,1]
# )
# df.snv.indel <- input$som
# genes_bed_path <- paste0(ROOT_DIR, '/data/gene_selection/genes.bed')
# genes.bed <- read.delim(gzfile(genes_bed_path), check.names=F)

mkMutProfileSnvIndel <- function(
   df.snv.indel, tumor.purity, mode,
   ignore.additional.evidence=F,
   rm.non.selected.genes=T, genes.bed=NULL,
   gene.identifier='ensembl_gene_id',
   keep.only.first.eff=T,
   verbose=T
){
   
   #--------- Sanity checks ---------#
   if(nrow(df.snv.indel)==0){ stop(return(NA)) }
   if(!(mode %in% c('germline','somatic'))){ stop("Please specify mode: 'germline', 'somatic'") }
   if(length(tumor.purity) != 1){ stop("tumor.purity must be a single numeric value") }


   #--------- Score annotations ---------#
   if(verbose){ message('Getting snpEff scores...') }
   # snpeff_score <- unlist(lapply(df.snv.indel$snpeff_eff, function(i){
   #    #i=snpeff.eff[[530]]
   #    if(grepl('&',i)){ i <- strsplit(i,'&')[[1]][1] }
   #    SCORING$snpeff[i]
   # }), use.names=F)
   snpeff_eff <- unlist(lapply(df.snv.indel$snpeff_eff, function(i){
      if(grepl('&',i)){ i <- strsplit(i,'&')[[1]][1] }
      return(i)
   }))
   snpeff_score <- unname(SCORING$snpeff)[ match(snpeff_eff,names(SCORING$snpeff)) ]
   snpeff_score[is.na(snpeff_score)] <- 0 ## If annotation not found in snpeff table, return 0

   if(verbose){ message('Getting ClinVar scores...') }
   # clinvar_score <- unlist(lapply(df.snv.indel$clinvar_sig, function(i){
   #    ifelse(is.na(i), 0, SCORING$clinvar[i])
   # }), use.names=F)
   clinvar_score <- unname(SCORING$clinvar)[ match(df.snv.indel$clinvar_sig,names(SCORING$clinvar)) ]
   clinvar_score[is.na(clinvar_score)] <- 0

   if(verbose){ message('Getting ENIGMA scores...') }
   # enigma_score <- unlist(lapply(df.snv.indel$enigma_sig, function(i){
   #    ifelse(is.na(i), 0, SCORING$enigma[i])
   # }), use.names=F)
   enigma_score <- unname(SCORING$enigma)[ match(df.snv.indel$enigma_sig,names(SCORING$enigma)) ]
   enigma_score[is.na(enigma_score)] <- 0

   sig_scores <- data.frame(
      snpeff_score, clinvar_score, enigma_score#,
      # cadd_int_score,
      # mcap_int_score,
      # scap_int_score
   )

   ## Calculate max score and which database it came from
   if(verbose){ message('Calculating max scores...') }
   db_names <- sapply(strsplit(colnames(sig_scores),'_'), `[`, 1, USE.NAMES = F)
   
   # max_sig_scores <- do.call(rbind, lapply(1:nrow(sig_scores), function(i){
   #    #i = 83
   #    r <- as.numeric(sig_scores[i,])
   #    max_score <- max(r, na.rm = T)
   #    
   #    if(max_score == 0){ max_score_origin <- 'none' }
   #    else{
   #       max_score_origin <- db_names[which(r == max_score)]
   #       max_score_origin <- paste(max_score_origin, collapse=',')
   #    }
   #    
   #    return(data.frame(max_score, max_score_origin))
   # }))
   
   max_sig_scores <- do.call(rbind, apply(sig_scores,1,function(i){
      max_score <- max(i)
      if(max_score == 0){ max_score_origin <- 'none' }
      else{
         max_score_origin <- db_names[which(i == max_score)]
         max_score_origin <- paste(max_score_origin, collapse=',')
      }
      return(data.frame(max_score, max_score_origin))
   }))
   
   sig_scores <- cbind(sig_scores, max_sig_scores)

   ### Modifiy out ###
   out <- cbind(
      df.snv.indel,
      sig_scores
   )

   #--------- Calculate adjusted tumor allele depth ---------#
   if(verbose){ message('Calculating adjusted tumor allele depths...') }
   if(ignore.additional.evidence){
      adj_tumor_ad_ref <- NA
      adj_tumor_ad_alt <- NA
      alt_exists <- NA
   } else {
      ## Calculate tumor purity adjusted AD
      ## Do outside of Map loop to take advantage of R's speed with vectors
      tumor_dp <- df.snv.indel$tumor_ad_ref + df.snv.indel$tumor_ad_alt

      normal_purity <- 1 - tumor.purity
      normal_dp <- tumor_dp * normal_purity
      normal_ad <- normal_dp/2 ## theoretical depth for each allele in normal tissue

      adj_tumor_ad_ref <- round(df.snv.indel$tumor_ad_ref - normal_ad)

      if(mode == 'germline'){
         adj_tumor_ad_alt <- round(df.snv.indel$tumor_ad_alt - normal_ad)

      } else if(mode == 'somatic'){
         ## Assumption: somatic variants called by strelka (somatic caller) come with the assumption
         ## that they belong exclusively to the tumor. Therefore, ALT_DP in somatic variants belongs
         ## exclusively to the tumor.
         adj_tumor_ad_alt <- df.snv.indel$tumor_ad_alt
      }

      alt_exists <- adj_tumor_ad_alt >= CUTOFFS$min.adj.tumor.ad.alt
   }

   ### Modifiy out ###
   out <- insColAfter(
      out,
      cbind(tumor_purity = tumor.purity, adj_tumor_ad_ref, adj_tumor_ad_alt, alt_exists),
      after = 'tumor_ad_alt'
   )

   #--------- Determine if ref is lost ---------#
   if(verbose){ message('Determining if ref is lost...') }
   if(ignore.additional.evidence){
      ad_diff_score <- NA
      ref_loss <- NA
   } else {
      ad_diff_score <- unlist(Map(function(adj_tumor_ad_ref, adj_tumor_ad_alt){
         ## ad_diff_score = alt/ref
         ## Taking into account where: ref = 0, ref < 0, and alt <= 0 (to prevent errors)
         if( is.na(adj_tumor_ad_ref) || is.na(adj_tumor_ad_alt) || adj_tumor_ad_alt <= 0 ){
            ad_diff_score <- 0
            
         } else if(adj_tumor_ad_ref < 0){
            ## alt/ref makes no sense when ref is negative.
            ## Set ref to 0: returns Inf or NaN
            ad_diff_score <- adj_tumor_ad_alt/0
            
         } else {
            ad_diff_score <- adj_tumor_ad_alt/adj_tumor_ad_ref
         }
         
         return(ad_diff_score)
      }, adj_tumor_ad_ref, adj_tumor_ad_alt, USE.NAMES=F))
      
      ref_loss <- ad_diff_score >= CUTOFFS$min.ad.diff.score
   }

   ### Modifiy out ###
   out <- insColAfter(
      out,
      cbind(ad_diff_score, ref_loss),
      after = 'alt_exists'
   )
   # }

   #--------- Get ENSG if if not exist ---------#
   ### Modifiy out ###
   if(gene.identifier=='snpeff_gene'){
      out$ensembl_gene_id <- geneNamesToEnsg(out$snpeff_gene)

      if(sum(is.na(out$ensembl_gene_id))!=0){
         message(sprintf(
            'Could not determine ENSEMBL gene id for %s/%s variants. Removing these...',
            sum(is.na(out$ensembl_gene_id)),
            nrow(out)
         ))
         out <- out[!is.na(out$ensembl_gene_id),]
      }
   }

   #--------- Get HGNC symbol ---------#
   ### Modifiy out ###
   out <- insColAfter(
      out,
      ensgToHgncSymbol(out$ensembl_gene_id),
      after = 'ensembl_gene_id',
      colname = 'hgnc_symbol'
   )

   #--------- Remove ENSG not provide by user (in genes bed file) ---------#
   ### Modifiy out ###
   if(rm.non.selected.genes){
      if(is.null(genes.bed)){
         warning('No genes.bed was not specified. Skipping removing non user selected genes')
      } else {
         if(verbose){ message('Removing non user selected genes...') }
         out <- out[out$ensembl_gene_id %in% genes.bed$ensembl_gene_id,]
      }
   }

   #--------- Deal with multiple sneff effs ---------#
   ### Modifiy out ###
   if(keep.only.first.eff){
      if(verbose){ message('Keeping only first snpeff_eff...') }
      out$snpeff_eff <- gsub('&.+$','',out$snpeff_eff)
   }

   #--------- Export ---------#
   return(out)
}
