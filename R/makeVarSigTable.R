library(VariantAnnotation)
options(stringsAsFactors=F)

base_dir <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'
if(dir.exists('/Users/lnguyen/')){
   base_dir <- paste0('/Users/lnguyen/', base_dir)
}

#========= Load databases =========#
readTableDefault <- function(path, ...){ read.table(path,sep='\t',header=T,...) }

var_sig_db_dir <- paste0(base_dir,'/scripts_main/hmfVariantAnnotation/data/variant_significance/')

clinvar_txt <- readTableDefault( paste0(var_sig_db_dir,'/clinvar_20181217.txt') )
enigma_txt <- readTableDefault( paste0(var_sig_db_dir,'/enigma_variants_20181221.txt') )

#========= Main =========#
makeVarSigTable <- function(vcf.file, mode){
   #vcf.file <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/variants_subset/CPCT02010419R_CPCT02010419T/CPCT02010419R_CPCT02010419T.multi.vcf.gz'
   
   vcf <- readVcf(vcf.file)
   rr <- rowRanges(vcf)
   
   ## Unselect rows with multiple ALT sequences
   if(!identical(which(lengths(rr$ALT) > 1), integer(0))){
      #if(verbose){ message('Some rows have multiple ALT sequences. These will be removed.') }
      sel_rows <- lengths(rr$ALT) == 1
   }
   
   vcf <- vcf[sel_rows]
   rr <- rr[sel_rows]
   
   ## Get coords
   df <- data.frame(
      chrom=as.character(seqnames(rr)),
      pos=start(rr),
      ref=as.character(rr$REF),
      alt=as.character(unlist(rr$ALT))
   )
   
   ## Get snpeff info from vcf
   ann <- info(vcf)$ANN
   ann <- lapply(ann, function(i){ str_split(i,'\\|')[[1]] })
   
   ann <- do.call(rbind, lapply(ann, function(i){
      data.frame(
         snpeff_eff = i[2],
         snpeff_gene = i[4],
         gene_id = i[5],
         transcript_id = i[7],
         hgvs_c = i[10]
      )
   }))
   
   ## Retrieve clinical significance
   getVarSigInDb <- function(db){
      #db <- enigma_txt
      
      do.call(rbind, lapply(1:nrow(df), function(i){
         #i=1
         r <- df[i,]
         sig <- db[
            db$chrom == r$chrom
            & db$pos == r$pos
            & db$ref == r$ref
            & db$alt == r$alt
            ,'sig']
         
         if(length(sig) > 0){
            return(sig)
         } else {
            return(NA)
         }
      }))
   }
   enigma_sig <- getVarSigInDb(enigma_txt)
   clinvar_sig <- getVarSigInDb(clinvar_txt)

   ## Extract genotype and allele depth of tumor sample
   extractFormat <- function(mode){
      
      if(!(mode %in% c('germline','somatic'))){ stop("Available modes are 'germline', 'somatic' ") }
      if(mode=='germline'){ 
         tumor_col <- 2
      } else if(mode=='somatic'){ 
         tumor_col <- 1
      }
      
      gt <- geno(vcf)$GT
      gt_tumor <- gt[,tumor_col]
      gt_tumor <- as.data.frame(do.call(rbind, str_split(gt_tumor,'/')))
      colnames(gt_tumor) <- c('gt1','gt2')
      
      ad <- as.data.frame(geno(vcf)$AD)
      ad_tumor <- unname(ad[,tumor_col])
      ad_tumor <-  as.data.frame(do.call(rbind, ad_tumor))
      colnames(ad_tumor) <- c('ad1','ad2')
      
      out <- cbind(gt_tumor, ad_tumor)
      return(out)
   }
   
   format <- extractFormat(mode)
   
   ## Output final table
   df <- cbind(df, ann, clinvar_sig, enigma_sig, format)
   
   return(df)
}

