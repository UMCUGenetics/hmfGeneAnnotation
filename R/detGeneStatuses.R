#!/usr/bin/env Rscript

#library(mltoolkit)
library(magrittr)
library(stringr)
options(stringsAsFactors = F)

#========= Paths =========#
base_dir <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'
if(dir.exists('/Users/lnguyen/')){
   base_dir <- paste0('/Users/lnguyen/', base_dir)
}

#========= Settings/constants =========#
HR_GENE_COORDS <- read.table(
   paste0(base_dir, 'scripts_main/hmfVariantAnnotation/data/gene_coords/hr_gene_coords_ensembl.bed'),
   sep='\t',header=F,col.names=c('chrom','start','end','entrez_id','gene','ensembl_id','aliases')
)

#--------- Cutoffs ---------#
CUTOFFS <- list(
   max.max.copy.number=0.2, ## full_gene_loss
   max.min.minor.allele.ploidy=0.2, ## loh
   min.cn.diff.in.gene=0.9, ## cn_break_in_gene
   min.adj.tumor.ad.alt=10, ## germ.alt_exists, som.alt_exists
   min.germ.ad.diff.score=1.5, ## germ.ref_loss
   hit.score.origin.min.score=6 ## hit_score_origin
)

#--------- Scoring ---------#
scoring_tables_dir <- paste0(base_dir,'/scripts_main/hmfVariantAnnotation/data/variant_significance/')
readScoringTables <- function(path){ read.table(path,sep='\t',header=T)[,c('ann','score')] }
SCORING <- list(
   ## databases
   snpeff=readScoringTables(paste0(scoring_tables_dir,'/snpeff_scoring.txt')),
   clinvar=readScoringTables(paste0(scoring_tables_dir,'/clinvar_scoring.txt')),
   enigma=readScoringTables(paste0(scoring_tables_dir,'/enigma_scoring.txt')),
   
   ## Use integer values for main evidance    
   full_gene_loss=30,
   loh=20,
   
   ## Additional evidence
   cn_break_in_gene_loh=0.2, ## When CN break happens with LOH, gene is more likely to be deficient
   cn_break_in_gene_germ_som=0.1,
   germ.ref_loss=0.1,
   
   ## Give bonus points when ALT AD is good
   germ.alt_exists=0.01,
   som.alt_exists=0.01
)

#--------- Var sig table headers ---------#
## Make headers for var_sig tables so that settings are recorded
SCORING_ss <- SCORING[!(names(SCORING) %in% c('snpeff','clinvar','enigma'))]
HEADERS <- list(
   cutoffs = sapply(1:length(CUTOFFS), function(i){ paste0(names(CUTOFFS[i]),'=',CUTOFFS[i]) }),
   scoring = sapply(1:length(SCORING_ss), function(i){ paste0(names(SCORING_ss)[i],'=',SCORING_ss[i]) })
)

#========= CNV =========#
calcCnvScores <- function(
   min.copy.number, max.copy.number, min.minor.allele.ploidy, 
   
   min.cn.diff.in.gene = CUTOFFS$min.cn.diff.in.gene, 
   max.max.copy.number = CUTOFFS$max.max.copy.number, 
   max.min.minor.allele.ploidy = CUTOFFS$max.min.minor.allele.ploidy, 
   
   show.raw=F
){
   # gene = cnv$Gene
   # min.copy.number = cnv$MinCopyNumber
   # max.copy.number = cnv$MaxCopyNumber
   # min.minor.allele.ploidy = cnv$MinMinorAllelePloidy
   
   ## Return empty variable if no CNVs present
   if(length(min.copy.number) == 0 
      & length(max.copy.number) == 0 & length(min.minor.allele.ploidy) == 0){
      cnv_scores <- NULL

   } else {
      cnv_scores <- Map(function(min.copy.number, max.copy.number, min.minor.allele.ploidy){
         ## Initiate scoring vector
         full_gene_loss=F
         loh=F
         cn_diff_in_gene=0
         cn_break_in_gene=F
         
         ## 'full gene loss' and 'loh' scores
         if(max.copy.number < max.max.copy.number){ full_gene_loss <- T }
         if(min.minor.allele.ploidy < max.min.minor.allele.ploidy){ loh <- T }
         
         cn_diff_in_gene <- abs(max.copy.number - min.copy.number) ## max should be > min, but do abs() just in case
         if(cn_diff_in_gene >= min.cn.diff.in.gene){ cn_break_in_gene <- T }
         
         return(data.frame(
            full_gene_loss, loh, cn_diff_in_gene, cn_break_in_gene
         ))
      }, min.copy.number, max.copy.number, min.minor.allele.ploidy, USE.NAMES=F)
      
      cnv_scores <- do.call(rbind,cnv_scores)
      
      if(show.raw){
         cnv_scores <- cbind(min.copy.number, max.copy.number, min.minor.allele.ploidy, cnv_scores)
      }
      
      #cnv_scores <- cbind(gene=gene,cnv_scores)
   }
   
   return(cnv_scores)
}

#========= SNV/indels =========#
calcSnvIndelScores <- function(snpeff.eff, clinvar.eff, enigma.eff, show.raw=F){
   # snpeff.eff=som$snpeff_eff
   # clinvar.eff=som$clinvar_sig
   # enigma.eff=som$enigma_sig
   # show.raw=T
   
   # snpeff.eff=germ$snpeff_eff
   # clinvar.eff=germ$clinvar_sig
   # enigma.eff=germ$enigma_sig
   # show.raw=T
   
   ## Return empty variable if SNVs/indels present
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
      snpeff_score <- unlist(lapply(snpeff.eff, function(i){
         #i=snpeff.eff[[986]]
         if(grepl('&',i)){
            i <- str_split(i,'&')[[1]]
            score <- max(unlist(lapply(i, function(j){
               SCORING$snpeff[SCORING$snpeff$ann == j,'score']
            })))
         } else {
            score <- SCORING$snpeff[SCORING$snpeff$ann == i,'score']
         }
         
         ## If annotation not found in snpeff table, return 0
         if(length(score) == 0 ){ score <- 0 }
         return(score)
      }))
      
      clinvar_score <- unlist(lapply(clinvar.eff, function(i){
         ifelse(
            is.na(i), 0, ## variants not found in db will have NA (from preprocessing with bash scripts)
            SCORING$clinvar[SCORING$clinvar$ann == i,'score']
         )
      }))
      
      enigma_score <- unlist(lapply(enigma.eff, function(i){
         ifelse(
            is.na(i), 0, 
            SCORING$enigma[SCORING$enigma$ann == i,'score']
         )
      }))
      
      out <- data.frame(snpeff_score, clinvar_score, enigma_score)
      
      ## Calculate max score and which database it came from 
      db_names <- str_remove_all(colnames(out),"_score")
      max_scores <- do.call(rbind, lapply(1:nrow(out), function(i){
         #i = 802
         r <- unlist(out[i,], use.names = F)
         max_score <- max(r)
         
         if(max_score == 0){ max_score_origin <- 'none' }
         else{
            max_score_origin <- db_names[which(r == max_score)]
            max_score_origin <- paste(max_score_origin, collapse=',')
         }
         
         return(data.frame(max_score, max_score_origin))
      }))
      
      out <- cbind(out, max_scores)
      
      if(show.raw){
         out <- cbind(snpeff.eff, clinvar.eff, enigma.eff, out)
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

#========= Format gene names =========#
getVarGene <- function(var.chrom, var.start, var.end=NULL, gene.coords=HR_GENE_COORDS){
   #var.chrom=1
   #var.start=28112808
   #var.start=28218035
   #var.end=28218036
   #var.end=28221257
   #var.end=28561768

   ## SNVs
   if(is.null(var.end)){
      isBetween <- function(x, l, r, incl.boundaries=T){
         if(incl.boundaries){ x >= l & x <= r }
         else { x > l & x < r }
      }

      gene <- gene.coords[
         var.chrom == gene.coords$chrom
         & isBetween(var.start, gene.coords$start, gene.coords$end)
         ,'gene']

      if(length(gene) > 1){
         ## greedy selection; in case an SNV occurs in a region where 2 genes overlap
         gene <- gene[1]
      }

      ## ranges; indels, SVs, CNVs
   } else {
      if(var.end < var.start){ stop('var.end must not be less than var.start') }

      nBasesOverlap <- function(start1, end1, start2, end2){
         min(end1, end2) - max(start1, start2)
      }

      gene.coords$var_n_bases_overlap <- unlist(Map(
         nBasesOverlap, var.start, var.end, gene.coords$start, gene.coords$end
      ))

      gene.coords <- gene.coords[gene.coords$var_n_bases_overlap > 0,]

      ## return gene immediately if only 1 match
      if(nrow(gene.coords_ss) == 1){
         gene <- gene.coords_ss$gene

         ## further processing if mulitple matches
      } else {
         ## which.max() is greedy; will select only the first match
         gene <- gene.coords[which.max(gene.coords$var_n_bases_overlap),'gene']
      }
   }

   if(length(gene) == 0){
      warning('No genes matched in gene.coords. Returning NA')
      gene <- NA
   }

   return(gene)
}

#========= Main =========#
detGeneStatuses <- function(
   sample.dir=paste0(sample.dir,'/'), ## Ensures that path is correct
   #sample.dir = paste0(base_dir,'/HMF_update/variants_subset/CPCT02010419R_CPCT02010419T/'),
   
   ## Default paths
   sample.name = str_split(basename(sample.dir), '_')[[1]][2],
   output.dir = paste0(sample.dir,'/gene_mut_profile/'), 
   cnv.table.path = paste0(sample.dir,'/',sample.name,'.purple.gene.cnv'),
   germ.table.path = paste0(sample.dir,'/germ_varsigs.txt'),
   som.table.path = paste0(sample.dir,'/som_varsigs.txt'),
   purity.table.path = paste0(sample.dir,'/',sample.name,'.purple.purity'),
   
   ## Options
   keep.only.first.eff=T,
   export.prelim.data=T,
   overwrite.mut.profiles=T,
   verbose=T
){
   ## Make output dir
   #output.dir <- paste0(sample.dir,'/gene_mut_profile/')
   if(!dir.exists(output.dir)){ dir.create(output.dir) }
   
   done_mark <- paste0(output.dir,'/done')
   if(file.exists(done_mark)){ file.remove(done_mark) }
   
   #--------- Load data ---------#
   if(verbose){ message('Loading data: ', sample.dir) }
   cnv <- read.table(cnv.table.path,sep='\t',header=T)
   germ <- read.table(germ.table.path,sep='\t',header=T)
   som <- read.table(som.table.path,sep='\t',header=T)
   tumor_purity <- read.table(purity.table.path, skip = 1)[,1]
   
   #::::::::: Testing :::::::::#
   # output.dir <- paste0(sample.dir,'/gene_mut_profile/')
   # if(!dir.exists(output.dir)){ dir.create(output.dir) }
   # 
   ## Load data
   # sample.dir = paste0(base_dir,'/HMF_update/variants_subset/CPCT02070360R_CPCT02070360T/')
   # sample.dir = paste0(base_dir,'/HMF_update/variants_subset/CPCT02020369R_CPCT02020369T/')
   # sample.dir = paste0(base_dir,'/HMF_update/variants_subset/CPCT02010815R_CPCT02010815T/')
   # sample.dir = paste0(base_dir,'/HMF_update/variants_subset/CPCT02020572R_CPCT02020572T/')
   # output.dir = paste0(sample.dir,'/gene_mut_profile/')
   # sample.name = str_split(basename(sample.dir), '_')[[1]][2]
   # 
   # sample.dir = paste0(base_dir,'/Ovarian_Organoids_Chris/vcf_subset/HGS_19_T27__HGS_19_B27/')
   # output.dir = paste0(sample.dir,'/gene_mut_profile/')
   # sample.name = 'HGS_19_T27__HGS_19_B27'
   # 
   # cnv <- read.table(paste0(sample.dir,'/',sample.name,'.purple.gene.cnv'),sep='\t',header=T)
   # germ <- read.table(paste0(sample.dir,'germ_varsigs.txt'),sep='\t',header=T)
   # som <- read.table(paste0(sample.dir,'som_varsigs.txt'),sep='\t',header=T)
   # tumor_purity <- read.table(paste0(sample.dir,'/',sample.name,'.purple.purity'), skip = 1)[,1]

   #--------- Assign numerical value to evidence ---------#
   mut_profile_paths <- list(
      cnv = paste0(output.dir,'/mut_profile_cnv.txt'),
      germ = paste0(output.dir,'/mut_profile_germ.txt'),
      som = paste0(output.dir,'/mut_profile_som.txt')
   )
   
   ## Calculate mut scores 
   if(!all(sapply(mut_profile_paths, file.exists)) | overwrite.mut.profiles){
      ## Fill in empty gene names
      if(verbose){ message('Filling in empty gene names...') }
      if(nrow(cnv) != 0){
         for(i in 1:nrow(cnv)){
            if(!is.na(cnv$Gene[i]) & nchar(cnv$Gene[i]) == 0){
               cnv$Gene[i] <- getVarGene(cnv$Chromosome[i], cnv$Start[i], cnv$End[i])
            }
         }
      }
      
      if(nrow(germ) != 0){
         for(i in 1:nrow(germ)){
            if(!is.na(germ$snpeff_gene[i]) & nchar(germ$snpeff_gene[i]) == 0){
               germ$snpeff_gene[i] <- getVarGene(germ$chrom[i], germ$pos[i])
            }
         }
      }
      
      if(nrow(som) != 0){
         for(i in 1:nrow(som)){
            if(!is.na(som$snpeff_gene[i]) & nchar(som$snpeff_gene[i]) != 0){
               som$snpeff_gene[i] <- getVarGene(som$chrom[i], som$pos[i])
            }
         }
      }
      
      #--------- Calculate variant scores ---------#
      ## Initiate main data structure
      gene_mut_profile <- list()
      
      ## Calculate scores
      if(verbose){ message('Calculating CNV scores...') }
      gene_mut_profile$cnv <- cbind(
         gene=cnv$Gene,
         calcCnvScores(
            cnv$MinCopyNumber, cnv$MaxCopyNumber, cnv$MinMinorAllelePloidy, 
            show.raw = T
         )
      )
      
      if(verbose){ message('Calculating germline/somatic mutation scores...') }
      ## Wrapper function for calcSnvIndelScores() and calcAdjTumorAd()
      makeGeneMutProfileGermSom <- function(df, mode){
         #df=som
         cbind(
            data.frame(gene = df$snpeff_gene),
            df[,c('chrom','pos','ref','alt','hgvs_c')],
            #df$chrom, df$pos, df$ref, df$alt, df$hgvs_c, df$gene_id,
            calcSnvIndelScores(df$snpeff_eff, df$clinvar_sig, df$enigma_sig, show.raw=T),
            calcAdjTumorAd(df$tumor_ad_ref, df$tumor_ad_alt, tumor_purity, mode = mode)
         )
      }
      
      gene_mut_profile$germ <- makeGeneMutProfileGermSom(germ, mode = 'germline')
      gene_mut_profile$som <- makeGeneMutProfileGermSom(som, mode = 'somatic')
      
      if(keep.only.first.eff){
         keepFirstEff <- function(df){
            #df=gene_mut_profile$germ
            df$snpeff.eff <- str_replace_all(df$snpeff.eff,'&.+$','')
            return(df)
         }
         
         gene_mut_profile$germ <- keepFirstEff(gene_mut_profile$germ)
         gene_mut_profile$som <- keepFirstEff(gene_mut_profile$som)
      }
      
      ## Export raw processed variant data before getting max effect
      ## This data is needed for identifying VUS's as pathogenic
      if(verbose){ message('Exporting raw mutation profile...') }
      if(export.prelim.data){
         write.tsv <- function(...){ write.table(..., sep='\t',quote=F,row.names=F) }
         
         write.tsv(gene_mut_profile$cnv, mut_profile_paths$cnv)
         write.tsv(gene_mut_profile$germ, mut_profile_paths$germ)
         write.tsv(gene_mut_profile$som, mut_profile_paths$som)
      }
      
   } else {
      gene_mut_profile <- lapply(mut_profile_paths, function(i){
         read.table(i,sep='\t',header=T)
      })
   }
   
   #--------- Make table containing 1 variant per type (i.e. cnv, germ, som) ---------#
   getGeneMaxEff <- function(df,greedy.max.eff=T){
      #df=gene_mut_profile$germ
      if(nrow(df) == 0){ out <- NA }
      else{
         genes <- unique(df$gene)
         #genes <- unique(df$gene_id)
         out <- do.call(rbind, lapply(genes, function(i){
            #i="BRE"
            
            df_ss <- df[ df$gene == i,]
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
   
   gene_statuses <- (function(){
      l <- list(
         cnv=gene_mut_profile$cnv,
         germ=getGeneMaxEff(gene_mut_profile$germ, greedy.max.eff=T),
         som=getGeneMaxEff(gene_mut_profile$som, greedy.max.eff=T)
      )
      
      union_genes <- unlist(lapply(l, function(i){
         if( is.data.frame(i) ){ i$gene }
         else { NA }
      }), use.names = F)
      
      union_genes <- union_genes[!is.na(union_genes) & nchar(union_genes) != 0] ## rm unknown genes and empty dataframes
      union_genes <- unique(sort(union_genes))
      
      ## Cols to select from gene_mut_profile
      cols <- list(
         cnv = c('min.minor.allele.ploidy','min.copy.number','max.copy.number','cn_diff_in_gene',
                 'full_gene_loss','loh','cn_break_in_gene'),
         
         germ = c('snpeff.eff','max_score','max_score_origin','adj_tumor_ad_ref',
                  'adj_tumor_ad_alt','alt_exists','ref_loss'),
         
         som  = c('snpeff.eff','max_score','max_score_origin','adj_tumor_ad_ref',
                  'adj_tumor_ad_alt','alt_exists')
      )
      
      out <- lapply(c('cnv','germ','som'), function(i){
         #i='cnv'
         #i='som'
         #i='germ'
         
         ## For df in gene_mut_profile with no variants, create df of NAs for errorless cbind
         if(!is.data.frame(l[[i]])){
            df <- data.frame(matrix( nrow=length(union_genes), ncol=length(cols[[i]]) ))
            colnames(df) <- cols[[i]]
            
         } else {
            df <- do.call(rbind, lapply(union_genes, function(j){
               #j=union_genes[79]
               df_temp <- l[[i]][ l[[i]]$gene == j, cols[[i]] ]
               
               ## If gene not in gene_mut_profile df, create df of NAs for errorless rbind (inner loop)
               if(nrow(df_temp) == 0){ df_temp[1,] <- rep(NA, length(cols[[i]])) }
               
               return(df_temp)
            }))
         }
         
         return(df)
      })
      
      names(out) <- c('','germ','som') ## Empty name for CNV results in no prefix with do.call(cbind, ...)
      out <- do.call(cbind, out)
      out <- cbind(gene = union_genes, out)
   })()
   
   #--------- Determine hit score, then deficient genes ---------#
   if(verbose){ message('Determining deficient genes...') }
   
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
      
      if(hit_score < CUTOFFS$hit.score.origin.min.score){
         hit_score_origin <- 'none'
      } else if(hit_score==SCORING$loh){
         hit_score_origin <- 'loh_only'
      } else {
         hit_score_origin <- names(hit_score_list)[hit_score_list==hit_score]
         hit_score_origin <- paste(hit_score_origin,collapse=',')
      }
      
      return(data.frame(is_def, def_type, hit_score, hit_score_origin))
   }))

   gene_statuses <- cbind(gene_statuses,def_status)
   
   ## Reorder columns for better readability
   sel_cols <- list(
      main=c("gene","is_def","def_type","hit_score","hit_score_origin"),
      
      cnv=c("min.minor.allele.ploidy","min.copy.number","max.copy.number",
            "cn_diff_in_gene","full_gene_loss","loh","cn_break_in_gene"),
      
      germ=c("germ.snpeff.eff","germ.max_score","germ.max_score_origin",
             "germ.alt_exists","germ.ref_loss","germ.adj_tumor_ad_ref","germ.adj_tumor_ad_alt"),
      
      som=c("som.snpeff.eff","som.max_score","som.max_score_origin",
            "som.alt_exists","som.adj_tumor_ad_ref","som.adj_tumor_ad_alt")
   )
   
   gene_statuses <- gene_statuses[,unlist(sel_cols)]
   
   ## Rename 
   colnames(gene_statuses) <- str_replace_all(colnames(gene_statuses),'[.]eff','_eff')
   
   #--------- Make expanded gene statuses table ---------#
   constant_cols <- unlist( sel_cols[c('main','cnv')] )
   variable_cols <- unique(colnames(gene_mut_profile$germ),colnames(gene_mut_profile$som))
   
   final_cols <- unname(c(constant_cols,variable_cols,'mut_origin'))
   
   gene_statuses_expanded <- lapply(1:nrow(gene_statuses), function(i){
      #i=9
      #print(i)
      r <- gene_statuses[i,]
      r[is.na(r)] <- 'NA' ## NA as a string allows this comparison: NA == NA; TRUE
      
      constant_rows <- r[, constant_cols ]
      
      if(r$hit_score_origin %in% c('full_gene_loss','loh_only','none')){
         out <- cbind(constant_rows, mut_origin=NA)
         
      } else {
         
         ## Get gene subset of gene_mut_profile$germ or $som
         getRelevantMuts <- function(type){
            #type='germ'
            gmp <- gene_mut_profile[[type]] %>% .[.$gene == r$gene,]
            
            gmp[is.na(gmp)] <- 'NA'
            
            if(nrow(gmp) == 0){ out <- data.frame() }
            else {
               if(type=='germ'){
                  out <- gmp %>% .[
                     .$max_score == r$germ.max_score
                     & .$ref_loss == r$germ.ref_loss
                     & .$alt_exists == r$germ.alt_exists
                     ,]
                  
               } else if(type=='som'){
                  out <- gmp %>% .[
                     .$max_score == r$som.max_score
                     & .$alt_exists == r$som.alt_exists
                     ,]
               }
               
               out$mut_origin <- type
            }

            return(out)
         }
         
         variable_rows_germ <- NULL
         variable_rows_som <- NULL
         
         if(grepl('germ',r$hit_score_origin)){ variable_rows_germ <- getRelevantMuts('germ') }
         if(grepl('som',r$hit_score_origin)){ variable_rows_som <- getRelevantMuts('som') }
         
         ## If only variable_rows_germ or variable_rows_som exists
         if( is.null(variable_rows_som) || nrow(variable_rows_som) == 0 ){
            variable_rows <- variable_rows_germ
            
         } else if( is.null(variable_rows_germ) || nrow(variable_rows_germ) == 0 ){
            variable_rows <- variable_rows_som
         
         ## If both exist   
         } else {
            missing_cols <- colnames(variable_rows_germ) %>% .[!(. %in% colnames(variable_rows_som))]
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
   gene_statuses_expanded$gene.1 <- NULL
   colnames(gene_statuses_expanded) <- str_replace_all(colnames(gene_statuses_expanded),'[.]eff','_eff')
   
   #--------- Export ---------#
   if(verbose){ message('Exporting final table...') }
   
   writeGeneStatuses <- function(df, out_path){
      ## Write header comments
      #out_path <- paste0(output.dir,'/',sample.name,'_gene_statuses.txt')
      if(file.exists(out_path)){ 
         file.remove(out_path)
         file.create(out_path)
      }
      
      writeHeaderLine <- function(text){
         text <- paste0('## ',text)
         write(text, out_path, append=T)
      }
      
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
   }
   
   writeGeneStatuses(gene_statuses,paste0(output.dir,'/',sample.name,'_gene_statuses.txt'))
   writeGeneStatuses(gene_statuses_expanded,paste0(output.dir,'/',sample.name,'_gene_statuses_expanded.txt'))
   
   ## Mark done
   file.create(done_mark)
   if(verbose){ message('Done') }
}

#========= Exec =========#
args <- commandArgs(trailingOnly=TRUE)
#detGeneStatuses(paste0(base_dir,'/HMF_update/variants_subset/CPCT02010419R_CPCT02010419T/'))
detGeneStatuses(args[1])
#detGeneStatuses(args[1], args[2])
