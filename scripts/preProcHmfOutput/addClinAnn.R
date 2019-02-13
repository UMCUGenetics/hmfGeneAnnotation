options(stringsAsFactors=F)

#========= Load data =========#
## Input file
args <- commandArgs(trailingOnly=TRUE)
variants_path <- args[1]
#variants_path <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020459R_CPCT02020459T/CPCT02020459R_CPCT02020459T.germ.txt.gz'

variants <- read.delim(gzfile(variants_path))

## Variant databases
# base_dir <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'
# if(dir.exists('/Users/lnguyen/')){
#    base_dir <- paste0('/Users/lnguyen/', base_dir)
# }
# var_sig_db_dir <- paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/variant_significance/')

ROOT_DIR <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/'
CLINVAR_DB_PATH <- paste0(ROOT_DIR,'/data/variant_significance/clinvar_20181217.txt')
clinvar_db <- read.delim(CLINVAR_DB_PATH)

ENIGMA_DB_PATH <- paste0(ROOT_DIR,'/data/variant_significance/enigma_variants_20181221.txt')
enigma_db <- read.delim(ENIGMA_DB_PATH)

## Remove row with multiple ALT
variants <- variants[grep(',',variants$alt,invert=T),]

## Get clinical significance
#df <- variants[variants$snpeff_gene == 'BRCA1',]
getVarSigInDb <- function(df, db){
   #db <- enigma_txt
   
   current_chrom <- 'none'
   do.call(rbind, lapply(1:nrow(df), function(i){
      #i=1
      r <- df[i,]
      
      if(r$chrom != current_chrom){
         current_chrom <<- r$chrom
         message('Retrieving clinical significance at chromosome: ',current_chrom)
      }
      
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

variants$clinvar_sig <- getVarSigInDb(variants, clinvar_db)
variants$enigma_sig <- getVarSigInDb(variants, enigma_db)

## Overwrite original file
message('Exporting new variants table...')
write.table(variants, gzfile(variants_path),row.names=F,quote=F,sep='\t')

################


# #========= Main =========#
# makeVarSigTable <- function(vcf.file, mode){
#    #vcf.file <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/variants_subset/CPCT02010419R_CPCT02010419T/CPCT02010419R_CPCT02010419T.multi.vcf.gz'
#    #vcf.file <- '/Users/lnguyen//hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/161022_HMFregCPCT_FR12244725_FR12244336_CPCT02010419/161022_HMFregCPCT_FR12244725_FR12244336_CPCT02010419.filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf.gz'
#    
#    vcf <- readVcf(vcf.file)
#    rr <- rowRanges(vcf)
#    
#    ## Unselect rows with multiple ALT sequences
#    if(!identical(which(lengths(rr$ALT) > 1), integer(0))){
#       #if(verbose){ message('Some rows have multiple ALT sequences. These will be removed.') }
#       sel_rows <- lengths(rr$ALT) == 1
#    }
#    
#    vcf <- vcf[sel_rows]
#    rr <- rr[sel_rows]
#    
#    ## Get coords
#    df <- data.frame(
#       chrom=as.character(seqnames(rr)),
#       pos=start(rr),
#       ref=as.character(rr$REF),
#       alt=as.character(unlist(rr$ALT))
#    )
#    
#    ## Get snpeff info from vcf
#    ann <- info(vcf)$ANN
#    ann <- lapply(ann, function(i){ str_split(i,'\\|')[[1]] })
#    
#    ann <- do.call(rbind, lapply(ann, function(i){
#       data.frame(
#          snpeff_eff = i[2],
#          snpeff_gene = i[4],
#          gene_id = i[5],
#          transcript_id = i[7],
#          hgvs_c = i[10]
#       )
#    }))
#    
#    ## Retrieve clinical significance
#    getVarSigInDb <- function(db){
#       #db <- enigma_txt
#       
#       do.call(rbind, lapply(1:nrow(df), function(i){
#          #i=1
#          r <- df[i,]
#          sig <- db[
#             db$chrom == r$chrom
#             & db$pos == r$pos
#             & db$ref == r$ref
#             & db$alt == r$alt
#             ,'sig']
#          
#          if(length(sig) > 0){
#             return(sig)
#          } else {
#             return(NA)
#          }
#       }))
#    }
#    enigma_sig <- getVarSigInDb(enigma_txt)
#    clinvar_sig <- getVarSigInDb(clinvar_txt)
# 
#    ## Extract genotype and allele depth of tumor sample
#    extractFormat <- function(mode){
#       
#       if(!(mode %in% c('germline','somatic'))){ stop("Available modes are 'germline', 'somatic' ") }
#       if(mode=='germline'){ 
#          tumor_col <- 2
#       } else if(mode=='somatic'){ 
#          tumor_col <- 1
#       }
#       
#       gt <- geno(vcf)$GT
#       gt_tumor <- gt[,tumor_col]
#       gt_tumor <- as.data.frame(do.call(rbind, str_split(gt_tumor,'/')))
#       colnames(gt_tumor) <- c('gt1','gt2')
#       
#       ad <- as.data.frame(geno(vcf)$AD)
#       ad_tumor <- unname(ad[,tumor_col])
#       ad_tumor <-  as.data.frame(do.call(rbind, ad_tumor))
#       colnames(ad_tumor) <- c('ad1','ad2')
#       
#       out <- cbind(gt_tumor, ad_tumor)
#       return(out)
#    }
#    
#    format <- extractFormat(mode)
#    
#    ## Output final table
#    df <- cbind(df, ann, clinvar_sig, enigma_sig, format)
#    
#    return(df)
# }
# 
