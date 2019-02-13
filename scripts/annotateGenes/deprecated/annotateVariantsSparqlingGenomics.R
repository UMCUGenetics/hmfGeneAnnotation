library(mltoolkit)
library(SPARQL)
library(stringr)
options(stringsAsFactors = F)
base_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/metadata/'

#========= Misc functions =========#
merge_by_rownames <- function(...){
   df <- merge(..., by = 'row.names')
   
   rownames(df) <- df[,1]
   df[,1] <- NULL
   
   return(df)
}

write.clipboard <- function(data){
   clip <- pipe("pbcopy", "w")                       
   write.table(data, file=clip)                               
   close(clip)
}

read.clipboard <- function(sep = '\t', header = T, ...){
   read.table(pipe("pbpaste"), sep = sep, header = header)
}

#========= Select patients =========#
hmf_response <- read.table(paste0(base_dir, '/brca_annotation/hmf_response.txt'))
selected_patients <- rownames(hmf_response)

#========= Annotate HR genes =========#
hmf_variants <- list(
   somatic = read.table(paste0(base_dir,'/hrd_gene_annotation/sparql_mc_hr_genes_somatic_snv_indel.txt')),
   cnv = read.table(paste0(base_dir,'/hrd_gene_annotation/sparql_mc_hr_genes_cnv.txt'))
)

## Rearrange list to root > patient > somatic,cnv
hmf_variants <- lapply(selected_patients, function(i){
   list(
      somatic = hmf_variants$somatic %>% .[.$sampleid == i,],
      cnv = hmf_variants$cnv %>% .[.$sampleid == i,]
   )
})
names(hmf_variants) <- selected_patients

#--------- Main ---------#
getVariantEffects <- function(patient){
   #patient = hmf_variants[['CPCT02070370T']]
   #patient = hmf_variants[['CPCT02030347T']]
   #patient = hmf_variants[['CPCT02170010T']]
   
   ## CNV
   if(nrow(patient$cnv) == 0){
      cnv_score <- data.frame(
         gene = character(), 
         full_gene_loss = integer(), 
         loh = integer(),
         cn_break_in_gene = integer()
      )
      
   } else {
      cnv_score <- as.data.frame(do.call(rbind, lapply(1:nrow(patient$cnv), function(i){
         #i = 1
         row <- patient$cnv[i,]
         
         ## Initiate scoring vector
         score <- data.frame(gene = row$gene, full_gene_loss = 0, loh = 0, cn_break_in_gene = 0)
         
         ## 'full gene loss' and 'loh' scores
         if(row$maxcopynumber < 0.1){ score['full_gene_loss'] <- 10 }
         if(row$minminoralleleploidy < 0.2){ score['loh'] <- 5 }
         
         ## cn_break_in_gene not used in final annotation, but will be handy to identify the cause of
         ## HRD as full gene loss for samples with only LOH
         if(round(row$mincopynumber - row$maxcopynumber, 0) != 0){ ## Use subtraction instead of division for speed
            score['cn_break_in_gene'] <- 1
         } 
         
         return(score)
      })))
   }
   
   ## somatic SNV, indel
   if(nrow(patient$somatic) == 0){ 
      snv_indel_score <- data.frame(gene = character(), coding_effect = integer())
   } else {
      snv_indel_score <- as.data.frame(do.call(rbind, lapply(1:nrow(patient$somatic), function(i){
         #i = 1
         row <- patient$somatic[i,]
         
         score <- data.frame(gene = row$gene, coding_effect = 0)
         
         if(row$worstcodingeffect == 'NONSENSE_OR_FRAMESHIFT'){
            score['coding_effect'] <- 5
         } else if (row$worstcodingeffect == 'MISSENSE'){
            score['coding_effect'] <- 3
         }
         
         return(score)
      })))
      
      ## Keep only the variant with the most pathogenic effect (highest score)
      snv_indel_score <- snv_indel_score[order(snv_indel_score$coding_effect, decreasing = T),]
      snv_indel_score <- snv_indel_score[!duplicated(snv_indel_score$gene),]
   }

   out <- merge(cnv_score, snv_indel_score, by = 'gene', all = T)
   out[is.na(out)] <- 0
   out$sum <- rowSums(subset(out, select=-c(gene,cn_break_in_gene)))
   
   return(out)
}

hmf_variant_effects <- lapply(1:length(hmf_variants), function(i){ 
   sample <- names(hmf_variants)[i]
   message('Determining variant effects for: ', i, ') ', names(hmf_variants)[i])
   variant_effects <- getVariantEffects(hmf_variants[[i]])
   cbind(sample = sample, variant_effects)
})
hmf_variant_effects <- do.call(rbind, hmf_variant_effects)






#--------- Check if annotations are the same as those from Arne ---------#
hmf_variant_effects_ss <- hmf_variant_effects[
   hmf_variant_effects$sum >= 9 & 
      (hmf_variant_effects$gene == 'BRCA1' | hmf_variant_effects$gene == 'BRCA2')
   ,]

brca_def_patients <- list(
   arne = rownames(hmf_response)[hmf_response$response != 'none'],
   sparql = hmf_variant_effects_ss$sample
)

brca_annotations <- do.call(rbind, lapply(unique(unlist(brca_def_patients)), function(i){
   data.frame(
      sample = i,
      arne = if(i %in% brca_def_patients$arne){ 1 } else { 0 },
      sparql = if(i %in% brca_def_patients$sparql){ 1 } else { 0 }
   )
}))

## Bind arne's annotation and sparql-based annotation
brca_annotations <- merge(
   brca_annotations, 
   cbind(sample=rownames(hmf_response), hmf_response), all.x=T
)
brca_annotations <- subset(brca_annotations, select=-c(SUM_BRCA1, SUM_BRCA2, SUM_new))


## Load arne's full annotations
arne_annotations <- read.table(paste0(base_dir,'/brca_annotation/prelim_output_HMF_DR10-update-fix_20180828.txt'), header=T, sep='\t')
arne_annotations <- subset(arne_annotations, select=c(
   Sample_ID_1, BRCA1_germ, BRCA2_germ, BRCA1_som, BRCA2_som,
   LOH_BRCA1, LOH_BRCA2, BRCA1_LOSS, BRCA2_LOSS, SUM, SUM_new
))
colnames(arne_annotations)[1] <- 'sample'

brca_annotations <- merge(brca_annotations,arne_annotations, all.x = T)
brca_annotations <- merge(brca_annotations,hmf_variant_effects_ss, all.x = T)

brca_annotations <- brca_annotations[order(brca_annotations$arne, decreasing = T),]

write.table(brca_annotations, '/Users/lnguyen/Desktop/test.txt')