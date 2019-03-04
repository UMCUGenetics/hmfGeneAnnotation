options(stringsAsFactors = F)

#========= Paths =========#
#--------- Package ---------#
ROOT_DIR='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/'
if(dir.exists('/Users/lnguyen/')){
   ROOT_DIR <- paste0('/Users/lnguyen/', ROOT_DIR)
}

## Paths, scoring, cutoffs, and options
source(paste0(ROOT_DIR,'/scripts/annotateGenes/detGeneStatuses_init.R'))

#--------- Input ---------#
args <- commandArgs(trailingOnly=TRUE)

out_dir <- args[1]
input_paths <- list(
   cnv = args[2],
   germ = args[3],
   som = args[4],
   purity = args[5]
)
sample_name <- args[6]

# sample_name='CPCT02050135R_CPCT02050135T'
# out_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/',sample_name)
# input_paths <- list(
#    cnv = paste0(out_dir,'/',sample_name,'.purple.gene.cnv'),
#    germ = paste0(out_dir,'/varsig/',sample_name,'_varsigs_germ.txt.gz'),
#    som = paste0(out_dir,'/varsig/',sample_name,'_varsigs_som.txt.gz'),
#    purity = paste0(out_dir,'/',sample_name,'.purple.purity')
# )

input <- list(
   cnv = read.delim(input_paths$cnv),
   germ = read.delim(input_paths$germ),
   som = read.delim(input_paths$som),
   purity = read.table(input_paths$purity, skip=1)[,1]
)

#========= Convert predictors to values =========#
#mut_profile_dir <- paste0(out_dir,'/mut_profile/')
#if(!dir.exists(mut_profile_dir)){ dir.create(mut_profile_dir) }

mut_profile_paths <- list(
   cnv = paste0(out_dir,'/mut_profile_cnv.txt.gz'),
   germ = paste0(out_dir,'/mut_profile_germ.txt.gz'),
   som = paste0(out_dir,'/mut_profile_som.txt.gz')
)

if(!all(sapply(mut_profile_paths, file.exists)) | OPTIONS$overwrite.mut.profiles){
   
   #--------- Core functions ---------#
   calcCnvScores <- function(
      min.copy.number, max.copy.number, min.minor.allele.ploidy, 
      
      min.cn.diff.in.gene = CUTOFFS$min.cn.diff.in.gene, 
      max.max.copy.number = CUTOFFS$max.max.copy.number, 
      min.min.copy.number = CUTOFFS$min.min.copy.number,
      max.min.minor.allele.ploidy = CUTOFFS$max.min.minor.allele.ploidy, 
      
      show.raw=F
   ){
      min.copy.number = input$cnv$min_copy_number
      max.copy.number = input$cnv$max_copy_number
      min.minor.allele.ploidy = input$cnv$min_minor_allele_ploidy
      
      ## Return empty variable if no CNVs present
      if(length(min.copy.number) == 0 
         & length(max.copy.number) == 0 & length(min.minor.allele.ploidy) == 0){
         cnv_scores <- NULL
         
      } else {
         cnv_scores <- Map(function(min.copy.number, max.copy.number, min.minor.allele.ploidy){
            ## Initiate scoring vector
            full_gene_loss=F
            loh=F
            amp=F
            cn_diff_in_gene=0
            cn_break_in_gene=F
            
            ## 'full gene loss' and 'loh' scores
            if(max.copy.number < max.max.copy.number){ full_gene_loss <- T }
            if(min.minor.allele.ploidy < max.min.minor.allele.ploidy){ loh <- T }
            
            cn_diff_in_gene <- abs(max.copy.number - min.copy.number) ## max should be > min, but do abs() just in case
            if(cn_diff_in_gene >= min.cn.diff.in.gene){ cn_break_in_gene <- T }
            
            if(min.copy.number >= min.min.copy.number){ amp <- T }
            
            return(data.frame(
               full_gene_loss, loh, cn_diff_in_gene, cn_break_in_gene, amp
            ))
         }, min.copy.number, max.copy.number, min.minor.allele.ploidy, USE.NAMES=F)
         
         cnv_scores <- do.call(rbind,cnv_scores)
         
         if(show.raw){
            cnv_scores <- cbind(
               min_copy_number = min.copy.number, 
               max_copy_number = max.copy.number, 
               min_minor_allele_ploidy = min.minor.allele.ploidy, 
               cnv_scores
            )
         }
      }
      
      return(cnv_scores)
   }
   
   calcSnvIndelScores <- function(snpeff.eff, clinvar.eff, enigma.eff, cadd.phred, cap.score, cap.type, show.raw=F){
      # type='germ'
      # snpeff.eff=input[[type]]$snpeff_eff
      # clinvar.eff=input[[type]]$clinvar_sig
      # enigma.eff=input[[type]]$enigma_sig
      # cadd.phred=input[[type]]$cadd_phred
      # cap.score=input[[type]]$cap_score
      # cap.type=input[[type]]$cap_type
      # show.raw=T
      
      ## Return empty variable if no SNVs/indels present
      if(
         (is.na(snpeff.eff) & is.na(clinvar.eff) & is.na(enigma.eff)) || 
         (length(snpeff.eff) == 0 & length(clinvar.eff) == 0 & length(enigma.eff) == 0)
      ){
         if(show.raw){
            out <- data.frame(
               snpeff.eff=vector(), clinvar.eff=vector(), enigma.eff=vector(),
               snpeff_score=vector(), clinvar_score=vector(), enigma_score=vector(),
               max_score=vector(), max_score_origin=vector()
            )
         } else {
            out <- data.frame(
               snpeff_score=vector(), clinvar_score=vector(), enigma_score=vector(),
               max_score=vector(), max_score_origin=vector()
            )
         }
         
      } else {
         ## Split into separate lapply loops for easy debugging
         if(OPTIONS$verbose){ message('  Getting snpEff scores...') }
         snpeff_score <- unlist(lapply(snpeff.eff, function(i){
            #i=snpeff.eff[[530]]
            if(grepl('&',i)){ i <- strsplit(i,'&')[[1]][1] }
            score <- SCORING$snpeff[i]
            
            ## If annotation not found in snpeff table, return 0
            if(is.na(score)){ score <- 0 }
            return(score)
         }), use.names=F)
         
         if(OPTIONS$verbose){ message('  Getting ClinVar scores...') }
         clinvar_score <- unlist(lapply(clinvar.eff, function(i){
            ifelse(is.na(i), 0, SCORING$clinvar[i])
         }), use.names=F)
         
         if(OPTIONS$verbose){ message('  Getting ENIGMA scores...') }
         enigma_score <- unlist(lapply(enigma.eff, function(i){
            ifelse(is.na(i), 0, SCORING$enigma[i])
         }), use.names=F)
         
         # if(OPTIONS$verbose){ message('  Calculating integer CADD scores...') }
         # cadd_int_score <- unlist(lapply(cadd.phred, function(i){
         #    ifelse(is.na(i) || i <= CUTOFFS$min.cadd.phred, 0, 5)
         # }), use.names=F)
         
         if(OPTIONS$verbose){ message('  Calculating integer MCAP scores...') }
         mcap_int_score <- unlist(Map(function(cap.score, cap.type){
            if(is.na(cap.type) || cap.type != 'mcap'){ 0 }
            else{ ifelse(cap.score >= CUTOFFS$min.mcap.score, 5, 0) }
         }, cap.score, cap.type), use.names=F)
         
         if(OPTIONS$verbose){ message('  Calculating integer SCAP scores...') }
         scap_int_score <- unlist(Map(function(cap.score, cap.type){
            #print(paste(cap.score, cap.type))
            if(is.na(cap.type) || cap.type == 'mcap'){ 0 }
            else { ifelse(cap.score >= CUTOFFS$scap.cutoffs[cap.type], 5, 0) }
         }, cap.score, cap.type), use.names=F)
         
         #out <- data.frame(snpeff_score, clinvar_score, enigma_score, cadd_int_score, mcap_int_score, scap_int_score)
         
         ## Exclude CADD score
         out <- data.frame(snpeff_score, clinvar_score, enigma_score, mcap_int_score, scap_int_score)
         
         ## Calculate max score and which database it came from
         db_names <- sapply(strsplit(colnames(out),'_'), `[`, 1, USE.NAMES = F)
         max_scores <- do.call(rbind, lapply(1:nrow(out), function(i){
            #i = 83
            r <- unlist(out[i,], use.names = F)
            max_score <- max(r, na.rm = T)
            
            if(max_score == 0){ max_score_origin <- 'none' }
            else{
               max_score_origin <- db_names[which(r == max_score)]
               max_score_origin <- paste(max_score_origin, collapse=',')
            }
            
            return(data.frame(max_score, max_score_origin))
         }))
         
         out <- cbind(out, max_scores)
         
         if(show.raw){
            out <- cbind(
               snpeff_eff=snpeff.eff, 
               clinvar_eff=clinvar.eff, 
               enigma_eff=enigma.eff,
               cadd_phred=cadd.phred,
               cap_score=cap.score,
               cap_type=cap.type,
               out
            )
         }
      }
      return(out)
   }
   
   calcAdjTumorAd <- function(tumor.ad.ref, tumor.ad.alt, tumor.purity, mode=NULL, show.raw=F){
      # tumor.ad.ref = germ$tumor_ad_ref
      # tumor.ad.alt = germ$tumor_ad_alt
      # 
      # tumor.ad.ref = som$tumor_ad_ref
      # tumor.ad.alt = som$tumor_ad_alt
      
      # tumor.purity = tumor_purity
      
      if(!(mode %in% c('germline','somatic'))){ stop("Please specify mode: 'germline', 'somatic'") }
      if(length(tumor.purity) != 1){ stop("tumor.purity must be a single numeric value") }
      
      ## Calculate tumor purity adjusted AD 
      ## Do outside of Map loop to take advantage of R's speed with vectors
      tumor_dp <- tumor.ad.ref + tumor.ad.alt
      
      normal_purity <- 1 - tumor.purity
      normal_dp <- tumor_dp * normal_purity
      normal_ad <- normal_dp/2 ## theoretical depth for each allele in normal tissue
      
      adj_tumor_ad_ref <- round(tumor.ad.ref - normal_ad)
      
      if(mode == 'germline'){
         adj_tumor_ad_alt <- round(tumor.ad.alt - normal_ad)
         
      } else if(mode == 'somatic'){
         ## Assumption: somatic variants called by strelka (somatic caller) come with the assumption 
         ## that they belong exclusively to the tumor. Therefore, ALT_DP in somatic variants belongs 
         ## exclusively to the tumor.
         adj_tumor_ad_alt <- tumor.ad.alt
      }
      
      alt_exists <- adj_tumor_ad_alt >= CUTOFFS$min.adj.tumor.ad.alt
      #alt_exists[is.na(alt_exists)] <- F
      
      if(show.raw){
         out <- data.frame(tumor.purity, tumor.ad.ref, tumor.ad.alt, 
                           adj_tumor_ad_ref, adj_tumor_ad_alt, alt_exists)
      } else {
         out <- data.frame(adj_tumor_ad_ref, adj_tumor_ad_alt, alt_exists)
      }
      
      ## Post-processing
      if(mode == 'germline'){
         ## determine if ref is lost in germline
         out$ad_diff_score <- unlist(Map(function(adj_tumor_ad_ref, adj_tumor_ad_alt){
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
         
         out$ref_loss <- out$ad_diff_score >= CUTOFFS$min.germ.ad.diff.score
      }
      
      return(out)
   }
   
   #--------- Exec ---------#
   ## Initiate main data structure
   mut_profile <- list()
   
   ## CNV
   if(OPTIONS$verbose){ message('Calculating CNV scores...') }
   mut_profile$cnv <- cbind(
      input$cnv[,c('ensembl_gene_id','hgnc_symbol')],
      calcCnvScores(
         input$cnv$min_copy_number, input$cnv$max_copy_number, input$cnv$min_minor_allele_ploidy, 
         show.raw = T
      )
   )
   
   ## Wrapper function for calcSnvIndelScores() and calcAdjTumorAd()
   makeGeneMutProfileGermSom <- function(df, mode){
      cbind(
         df[,c('ensembl_gene_id','snpeff_gene','chrom','pos','ref','alt','hgvs_c')],
         calcSnvIndelScores(df$snpeff_eff, df$clinvar_sig, df$enigma_sig, df$cadd_phred, df$cap_score, df$cap_type, show.raw=T),
         calcAdjTumorAd(df$tumor_ad_ref, df$tumor_ad_alt, input$purity, mode = mode)
      )
   }
   
   if(OPTIONS$verbose){ message('Calculating germline mutation scores...') }
   mut_profile$germ <- makeGeneMutProfileGermSom(input$germ, mode = 'germline')
   
   if(OPTIONS$verbose){ message('Calculating somatic mutation scores...') }
   mut_profile$som <- makeGeneMutProfileGermSom(input$som, mode = 'somatic')
   
   colnames(mut_profile$germ)[1] <- 'ensembl_gene_id'
   colnames(mut_profile$som)[1] <- 'ensembl_gene_id'
   
   if(OPTIONS$keep.only.first.eff){
      keepFirstEff <- function(df){
         #df=mut_profile$germ
         df$snpeff_eff <- gsub('&.+$','',df$snpeff_eff)
         return(df)
      }
      
      mut_profile$germ <- keepFirstEff(mut_profile$germ)
      mut_profile$som <- keepFirstEff(mut_profile$som)
   }
   
   ## Export raw processed variant data before getting max effect
   ## This data is needed for identifying VUS's as pathogenic
   if(OPTIONS$verbose){ message('Exporting mutation profile...') }

   write.tsv <- function(...){ write.table(..., sep='\t',quote=F,row.names=F) }
   
   write.tsv(mut_profile$cnv, gzfile(mut_profile_paths$cnv))
   write.tsv(mut_profile$germ, gzfile(mut_profile_paths$germ))
   write.tsv(mut_profile$som, gzfile(mut_profile_paths$som))
   
} else {
   if(OPTIONS$verbose){ message('Mutation profiles exist. Loading tables...') }
   mut_profile <- lapply(mut_profile_paths, function(i){ read.delim(gzfile(i)) })
}




#========= Make table containing most pathogenic variant combination per gene =========#
if(OPTIONS$verbose){ message('Determining most pathogenic variant combination per gene...') }

getGeneMaxEff <- function(df, greedy.max.eff=T){
   #df=gene_mut_profile$germ
   if(nrow(df) == 0){ out <- NA }
   else{
      gene_ids <- unique(df$ensembl_gene_id)
      out <- do.call(rbind, lapply(gene_ids, function(i){
         #i="BRE"
         
         df_ss <- df[ df$ensembl_gene_id == i,]
         df_ss[which.max(df_ss$max_score),]
         
         if(greedy.max.eff){
            out <- df_ss[which.max(df_ss$max_score),]
         } else {
            max_max_score <- max(df_ss$max_score)
            out <- df_ss[df_ss$max_score == max_max_score, ]
         }
      }))
   }
   return(out)
}

## Load gene info to subset for selected genes
human_genes_enst2ensg <- read.delim(gzfile(GENES_ENST2ENSG), check.names=F)
genes_bed <- read.delim(GENES_BED, check.names=F)

## Cols to select from gene_mut_profile
sel_cols <- list(
   cnv = c('min_minor_allele_ploidy','min_copy_number','max_copy_number','cn_diff_in_gene',
           'full_gene_loss','loh','cn_break_in_gene','amp'),
   
   germ = c('snpeff_eff','clinvar_eff','enigma_eff','max_score','max_score_origin',
            'cadd_phred','cap_score','cap_type',
            'adj_tumor_ad_ref','adj_tumor_ad_alt','alt_exists','ref_loss'),
   
   som  = c('snpeff_eff','clinvar_eff','enigma_eff','max_score','max_score_origin',
            'cadd_phred','cap_score','cap_type',
            'adj_tumor_ad_ref','adj_tumor_ad_alt','alt_exists')
)

## Main
gene_statuses <- (function(){
   l <- list(
      cnv=mut_profile$cnv,
      germ=getGeneMaxEff(mut_profile$germ, greedy.max.eff=T),
      som=getGeneMaxEff(mut_profile$som, greedy.max.eff=T)
   )
   
   union_gene_ids <- unlist(lapply(l, function(i){
      if( is.data.frame(i) ){ i$ensembl_gene_id }
      else { NA }
   }), use.names = F)
   
   union_gene_ids <- union_gene_ids[!is.na(union_gene_ids) & nchar(union_gene_ids) != 0] ## rm unknown genes and empty dataframes
   union_gene_ids <- unique(union_gene_ids)
   
   out <- lapply(c('cnv','germ','som'), function(i){
      #i='cnv'
      #i='som'
      #i='germ'
      
      ## For df in gene_mut_profile with no variants, create df of NAs for errorless cbind
      if(!is.data.frame(l[[i]])){
         df <- data.frame(matrix( nrow=length(union_gene_ids), ncol=length(sel_cols[[i]]) ))
         colnames(df) <- sel_cols[[i]]
         
      } else {
         df <- do.call(rbind, lapply(union_gene_ids, function(j){
            #j=union_gene_ids[79]
            df_temp <- l[[i]][ l[[i]]$ensembl_gene_id == j, sel_cols[[i]] ]
            
            ## If gene id not in mut_profile df, create df of NAs for errorless rbind (inner loop)
            if(nrow(df_temp) == 0){ df_temp[1,] <- rep(NA, length(sel_cols[[i]])) }
            
            return(df_temp)
         }))
      }
      
      return(df)
   })
   
   names(out) <- c('','germ','som') ## Empty name for CNV results in no prefix with do.call(cbind, ...)
   out <- do.call(cbind, out)
   
   ## Add gene names and sort
   out <- cbind(
      ensembl_gene_id = union_gene_ids, 
      hgnc_symbol = human_genes_enst2ensg[match(union_gene_ids, human_genes_enst2ensg$ensembl_gene_id),'hgnc_symbol'],
      out)
   
   out <- out[order(out$hgnc_symbol),]
   
   return(out)
})()

## Remove genes not in gene selection
gene_statuses <- gene_statuses[gene_statuses$ensembl_gene_id %in% genes_bed$ensembl_gene_id,]


#========= Determine hit score, then deficient genes =========#
if(OPTIONS$verbose){ message('Determining deficient genes...') }

## IS_DEF_MIN_HIT_SCORE is ordered by descending likelihood of causing gene deficiency
IS_DEF_MIN_HIT_SCORE <- c(
   'full_gene_loss' = SCORING$full_gene_loss,
   
   ## loh & germ.max_score==5 & germ.alt_exists & (germ.ref_loss | cn_break_in_gene)
   'loh+germ' = 
      SCORING$loh + 5 + 
      SCORING$germ.alt_exists + min(SCORING$germ.ref_loss, SCORING$cn_break_in_gene_loh),
   
   ## loh & som.max_score==5 & som.alt_exists
   'loh+som' = SCORING$loh + 5 + SCORING$som.alt_exists,
   
   ## germ.max_score==5 & som.max_score==5 & germ.alt_exists & (germ.ref_loss | cn_break_in_gene) & som.alt_exists )
   'germ+som' = 5 + 5 + 
      SCORING$germ.alt_exists + min(SCORING$germ.ref_loss, SCORING$cn_break_in_gene_loh) + 
      SCORING$som.alt_exists
)

calcHitScores <- function(
   full_gene_loss=0, loh=0, germ.max_score=0, som.max_score=0,
   cn_break_in_gene=0, germ.ref_loss=0, germ.alt_exists=0, som.alt_exists=0, 
   output=NULL
){
   out <- c(
      'full_gene_loss' = ifelse(full_gene_loss, SCORING$full_gene_loss, 0),
      
      'loh+germ' =
         ifelse(loh, SCORING$loh, 0) +
         germ.max_score +
         ifelse(germ.alt_exists, SCORING$germ.alt_exists, 0) +
         
         ## *** Ensures that germ and som mutations have equal weighting, since 'som.ref_loss' does
         ## not exist/is not calculated (as it is unreliable), while 'germ.ref_loss' does exist
         max(
            ifelse(cn_break_in_gene, SCORING$cn_break_in_gene_loh, 0),
            ifelse(germ.ref_loss, SCORING$germ.ref_loss, 0)
         ),
      
      'loh+som' =
         ifelse(loh, SCORING$loh, 0) +
         som.max_score +
         ifelse(som.alt_exists, SCORING$som.alt_exists, 0) +
         ifelse(cn_break_in_gene, SCORING$cn_break_in_gene_loh, 0),
      
      'germ+som' =
         germ.max_score +
         ifelse(germ.alt_exists, SCORING$germ.alt_exists, 0) +
         max(
            ifelse(cn_break_in_gene, SCORING$cn_break_in_gene_loh, 0),
            ifelse(germ.ref_loss, SCORING$germ.ref_loss, 0)
         ) +
         som.max_score +
         ifelse(som.alt_exists, SCORING$som.alt_exists, 0)
   )
   
   if(is.null(output)){ 
      return(out) 
   } else {
      return(out[output])
   }
}

def_status <- do.call(rbind, lapply(1:nrow(gene_statuses), function(i){
   #i=3
   r <- gene_statuses[i,]
   r[is.na(r)] <- 0 ## treat NA's as 0 to avoid errors
   
   hit_score_list <- calcHitScores( 
      r$full_gene_loss, r$loh, r$germ.max_score, r$som.max_score,
      r$cn_break_in_gene, r$germ.ref_loss, r$germ.alt_exists, r$som.alt_exists
   )
   
   ## Determine is_def
   def_type <- 'none'
   is_def <- FALSE
   for(i in names(hit_score_list)){
      if(hit_score_list[[i]] >= IS_DEF_MIN_HIT_SCORE[[i]]){
         def_type <- i
         is_def <- TRUE
         break
      }
   }
   
   ## Return hit_score
   ## For determining likely cause of HRD when ClinVar/ENIGMA/SnpEff score is <= 4
   hit_score <- max(hit_score_list)
   
   if(hit_score < CUTOFFS$min.hit.score){
      hit_type <- 'none'
   } else if(hit_score==SCORING$loh){
      hit_type <- 'loh_only'
   } else {
      hit_type <- names(hit_score_list)[hit_score_list==hit_score]
      hit_type <- paste(hit_type,collapse=',')
   }
   
   return(data.frame(is_def, def_type, hit_score, hit_type))
}))

gene_statuses <- cbind(gene_statuses,def_status)

## Reorder columns for better readability
sel_cols <- (function(){ 
   l <- list(main=c('hgnc_symbol','ensembl_gene_id','is_def','def_type','hit_score','hit_type'))
   l <- c(l, sel_cols)
   return(l)
})()

sel_cols$germ <- paste0('germ.',sel_cols$germ)
sel_cols$som <- paste0('som.',sel_cols$som)

gene_statuses <- gene_statuses[,unlist(sel_cols)]

#========= Make expanded gene statuses table =========#
if(OPTIONS$verbose){ message('Retrieving variant combinations with maximum hit score per gene...') }

constant_cols <- unlist( sel_cols[c('main','cnv')], use.names=F )
variable_cols <- unique(colnames(mut_profile$germ),colnames(mut_profile$som))

final_cols <- unname(c(constant_cols,variable_cols,'mut_origin'))

gene_statuses_expanded <- lapply(1:nrow(gene_statuses), function(i){
   #i=297
   #print(i)
   r <- gene_statuses[i,]
   r[is.na(r)] <- 'NA' ## NA as a string allows this comparison: NA == NA; TRUE
   
   constant_rows <- r[, constant_cols ]
   
   if(r$hit_type %in% c('full_gene_loss','loh_only','none')){
      out <- cbind(constant_rows, mut_origin=NA)
      
   } else {
      
      ## Get gene subset of mut_profile$germ or $som
      getRelevantMuts <- function(type){
         #type='germ'
         gmp <- mut_profile[[type]][ mut_profile[[type]]$ensembl_gene_id == r$ensembl_gene_id, ]
         
         gmp[is.na(gmp)] <- 'NA'
         
         if(nrow(gmp) == 0){ out <- data.frame() }
         else {
            if(type=='germ'){
               out <- gmp[
                  gmp$max_score == r$germ.max_score
                  & gmp$ref_loss == r$germ.ref_loss
                  & gmp$alt_exists == r$germ.alt_exists
                  ,]
               
            } else if(type=='som'){
               out <- gmp[
                  gmp$max_score == r$som.max_score
                  & gmp$alt_exists == r$som.alt_exists
                  ,]
            }
            
            out$mut_origin <- type
         }
         
         return(out)
      }
      
      variable_rows_germ <- NULL
      variable_rows_som <- NULL
      
      if(grepl('germ',r$hit_type)){ variable_rows_germ <- getRelevantMuts('germ') }
      if(grepl('som',r$hit_type)){ variable_rows_som <- getRelevantMuts('som') }
      
      ## If only variable_rows_germ or variable_rows_som exists
      if( is.null(variable_rows_som) || nrow(variable_rows_som) == 0 ){
         variable_rows <- variable_rows_germ
         
      } else if( is.null(variable_rows_germ) || nrow(variable_rows_germ) == 0 ){
         variable_rows <- variable_rows_som
         
         ## If both exist   
      } else {
         missing_cols <- colnames(variable_rows_germ)[!(colnames(variable_rows_germ) %in% colnames(variable_rows_som))]
         variable_rows_som[,missing_cols] <- NA
         variable_rows_som <- variable_rows_som[,colnames(variable_rows_germ)]
         
         variable_rows <- rbind(variable_rows_germ,variable_rows_som)
      }
      
      constant_rows <- constant_rows[rep(1,nrow(variable_rows)),]
      
      out <- cbind(constant_rows,variable_rows)
   }
   
   if(nrow(out) != 0){
      missing_final_cols <- final_cols[!(final_cols %in% colnames(out))]
      out[,missing_final_cols] <- 'NA'
      out <- out[,final_cols]
   } else { ## For when no germ/som variants exist
      out <- NULL
   }
   
   return(out)
})

gene_statuses_expanded <- do.call(rbind,gene_statuses_expanded)
gene_statuses_expanded$ensembl_gene_id.1 <- NULL




#========= Export =========#
if(OPTIONS$verbose){ message('Exporting final table...') }

gene_statuses_paths <- list(
   main = paste0(out_dir,'/gene_statuses.txt'),
   expanded = paste0(out_dir,'/gene_statuses_expanded.txt')
)

SCORING_ss <- SCORING[!(names(SCORING) %in% c('snpeff','clinvar','enigma'))]
HEADERS <- list(
   cutoffs = sapply(1:length(CUTOFFS), function(i){ paste0(names(CUTOFFS[i]),'=',CUTOFFS[i]) }),
   scoring = sapply(1:length(SCORING_ss), function(i){ paste0(names(SCORING_ss)[i],'=',SCORING_ss[i]) })
)

writeGeneStatuses <- function(df, out_path){
   ## Write header comments
   #out_path <- paste0(out_dir,'/gene_statuses.txt')
   if(file.exists(out_path)){ file.remove(out_path) }
   file.create(out_path)
   
   writeHeaderLine <- function(text){
      text <- paste0('## ',text)
      write(text, out_path, append=T)
   }
   
   writeHeaderLine(paste0("SAMPLE_NAME: ",sample_name))
   writeHeaderLine("")
   
   writeHeaderLine("CUTOFFS")
   for(i in HEADERS$cutoffs){ writeHeaderLine(i) }
   writeHeaderLine("")
   
   writeHeaderLine("SCORING")
   for(i in HEADERS$scoring){ writeHeaderLine(i) }
   writeHeaderLine("")
   
   ## Export final table
   suppressWarnings(  ## Prevent 'appending column names to file' warning
      write.table(df, out_path, sep='\t',quote=F,row.names=F,append=T)
   )
   
   ## Compress final table
   system(sprintf("gzip -c %s > %s.gz && rm %s", out_path, out_path, out_path))
}

writeGeneStatuses(gene_statuses, gene_statuses_paths$main)
writeGeneStatuses(gene_statuses_expanded, gene_statuses_paths$expanded)


