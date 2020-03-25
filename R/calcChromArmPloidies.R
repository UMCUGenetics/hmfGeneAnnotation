#' Calculate overall chrom arm copy numbers
#'
#' @description This function first rounds copy numbers (CN) to integers so that CN segments can 
#' be grouped together. Per chrom arm, the coverage of each CN category is calculated (i.e. 
#' cumulative segment size). The chrom arm CN is (roughly) defined as the CN category with the
#' highest cumulative segment size
#'
#' @param purple.cnv.file Path to purple cnv file
#' @param out.file Path to output file. If NULL, returns a named vector
#' @param min.rel.cum.segment.size If a chrom arm has a CN category that covers >0.5 (i.e 50%; default) 
#' of a chrom arm, this CN is the copy number of the arm
#' @param max.rel.cum.segment.size.diff This value (default 0.1) determines whether which CN
#' categories are considered to cover equal lengths of the chrom arm. For example, (by default) 2
#' CN categories covering 0.40 and 0.31 of a chrom arm are considered equally contributing. When 
#' these CN categories have similar cumulative segment size as the one with the highest, if one of 
#' these have the same CN as the genome CN, return the genome CN. Otherwise, simply return the one 
#' with the highest segment support (as is done above).
#' @param chrom.arm.names A character vector in the form c('1p','1q','2p','2q', ...). The default 
#' 'auto' means that the human chromosome arm names are used. Note that chroms 13, 14, 15, 21, 22
#' are considered to only have the long (i.e. q) arm.
#' @param chrom.arm.split.method Which method to determine the chromosome arm coords? If 'hmf', uses
#' 'method' column from purple cnv file to determine centromere positions (i.e. p/q arm split point).
#' If 'gap', uses the a (processed) gap.txt.gz table from the UCSC genome browser to determine
#' centromere positions. These 2 methods should in theory be identical, unless the HMF pipeline code
#' changes.
#' @param verbose Show progress messages?
#'
#' @return A named vector of chrom arm copy numbers, or specified writes a table to out.file if 
#' specified
#' @export
#'
#' @examples
calcChromArmPloidies <- function(
   purple.cnv.file, out.file=NULL, 
   min.rel.cum.segment.size=0.5, max.rel.cum.segment.size.diff=0.1,
   chrom.arm.split.method='hmf', 
   centromere.positions.path=CENTROMERE_POSITIONS, one.armed.chroms=ONE_ARMED_CHROMS,
   chrom.arm.names='auto', 
   verbose=T
){
   
   #purple.cnv.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/data//CPCT02210082T.purple.cnv'
   #purple.cnv.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/data/CPCT02410017T.purple.cnv'
   #purple.cnv.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-DR047/data//160704_HMFregCPCT_FR12244557_FR12244595_CPCT02110002/CPCT02110002T.purple.cnv'
   
   cnv <- read.delim(purple.cnv.file, check.names=F, stringsAsFactors=F)
   colnames(cnv) <- sub('#','',colnames(cnv))
   
   #--------- Pre-calculations ---------#
   if(chrom.arm.split.method=='hmf'){
      cnv <- cnv[,c('chromosome','start','end','copyNumber','segmentStartSupport','segmentEndSupport','method')]
   }
   colnames(cnv)[1:4] <- c('chrom','start','end','copy_number')
   cnv$chrom <- gsub('chr','',cnv$chrom)
   
   cnv$segment_size <- (cnv$end - cnv$start) + 1
   cnv$copy_number[cnv$copy_number < 0] <- 0 ## Ceiling negative values
   cnv$copy_number_int <- round(cnv$copy_number)
   
   #--------- Split by chromosome and arm ---------#
   if(verbose){ message('Splitting CN segments by chromosome...') }
   cnv_split_pre <- split(cnv, cnv$chrom)
   cnv_split_pre <- cnv_split_pre[unique(cnv$chrom)] ## maintain chrom order
   
   if(verbose){ message('Splitting CN segments by arm...') }
   if(chrom.arm.split.method=='hmf'){
      cnv_split <- lapply(cnv_split_pre,function(i){
         #i <- cnv_split_pre[[22]]
         if(!('LONG_ARM' %in% i$method)){ ## One armed chromosomes are marked by LONG_ARM in method
            p_arm_end_index <- which(i$segmentEndSupport=='CENTROMERE')
            l <- list(
               p=i[1:p_arm_end_index,],
               q=i[(p_arm_end_index+1):nrow(i),]
            )
         } else {
            ## Consider one armed chroms as only having q arm
            l <- list(
               q=i
            )
         }
         names(l) <- paste0(i[1,'chrom'],names(l)) ## Make list names in the form 1p,1q,etc before flattening list
         return(l)
      })
   }
   
   if(chrom.arm.split.method=='gap'){
      centro_pos <- read.delim(centromere.positions.path, stringsAsFactors=F)
      centro_pos <- structure(centro_pos$pos, names=centro_pos$chrom)
      
      cnv_split <- lapply(cnv_split_pre,function(i){
         #i <- cnv_split_pre[['11']]
         chrom <- i[1,'chrom']
         if(!(chrom %in% one.armed.chroms)){
            q_arm_start_index <- which.min(i$end < centro_pos[chrom])
            
            l <- list(
               p=i[1:(q_arm_start_index-1),],
               q=i[q_arm_start_index:nrow(i),]
            )
         } else {
            l <- list(
               q=i
            )
         }
         names(l) <- paste0(i[1,'chrom'],names(l)) ## Make list names in the form 1p,1q,etc before flattening list
         return(l)
      })
   }
   
   ## Flatten 2-level list into 1-level list
   names(cnv_split) <- NULL
   cnv_split <- do.call(c, cnv_split)
   rm(cnv_split_pre)
   
   #--------- Calc arm ploidies ---------#
   if(verbose){ message('Calculating preliminary arm ploidies...') }
   cn_segment_support <- lapply(cnv_split, function(i){
      #i=cnv_split$`17q`
      
      df <- aggregate(i$segment_size, by=list(i$copy_number_int), FUN=sum)
      colnames(df) <- c('copy_number_int','cum_segment_size')
      
      df$cum_segment_size_rel <- df$cum_segment_size / sum(df$cum_segment_size)
      
      #df[which.max(df$cum_segment_size_rel),]
      df[order(df$cum_segment_size_rel, decreasing=T),]
   })
   
   ## CN with most frequent total segment support is preliminary CN
   arm_cn_prelim <- unlist(lapply(cn_segment_support, function(i){ i[1,'copy_number_int'] }))
   
   genome_cn <- as.numeric(names(
      sort(table(arm_cn_prelim),decreasing=T)
   )[1])
   
   if(verbose){ message('Calculating final arm ploidies...') }
   arm_cn <- unlist(lapply(cn_segment_support, function(i){
      #i=cn_segment_support[['12p']]
      
      ## E.g. >=50% of arm has CN of 2, then this is the CN
      if( i[1,'cum_segment_size_rel'] >= min.rel.cum.segment.size ){
         return(i[1,'copy_number_int'])
      }
      
      #' When multiple CNs have similar segment support as the one with the highest, if one 
      #' of these have the same CN as the genome CN, return the genome CN. Otherwise, simply return 
      #' the one with the highest segment support (as is done above)
      i$diffs <- i[1,'cum_segment_size_rel'] - i[,'cum_segment_size_rel']
      cn_doubt <- i[i$diffs < max.rel.cum.segment.size.diff,'copy_number_int']
      
      if(any(cn_doubt==genome_cn)){
         return(genome_cn)
      } else {
         return(i[1,'copy_number_int'])
      }
   }))
   
   ploidies <- c(arm_cn, genome=genome_cn)
   
   #--------- Ensure consistent output (e.g. for when CN data is missing) ---------#
   if(verbose){ message('Preparing output...') }
   if(is.null(chrom.arm.names)){
      return(ploidies)
   }
   
   if(chrom.arm.names=='auto'){
      chrom_arm_names <- paste0(rep(c(1:22,'X'),each=2),c('p','q'))
      chrom_arm_names <- chrom_arm_names[!(chrom_arm_names %in% paste0(one.armed.chroms,'p'))] ## Keep only q arm for one arm chromosomes
      chrom_arm_names <- c(chrom_arm_names,'genome')
   } else {
      chrom_arm_names <- c(chrom.arm.names,'genome')
   }
   
   out <- structure(rep(NA,length(chrom_arm_names)), names=chrom_arm_names)
   ploidies <- ploidies[names(ploidies) %in% chrom_arm_names] ## Rm unwanted chroms/arms (e.g. mitochondria)
   out[names(ploidies)] <- ploidies
   
   if(is.null(out.file)){ 
      return(out) 
   } else {
      out <- data.frame(chrom=names(out),ploidy=out,row.names=NULL)
      if(verbose){ message('Writing output...') }
      write.tsv(out, out.file)
   }
   
}

# out <- calcChromArmPloidies('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/ploidy_estimator/data/CPCT02410017T.purple.cnv')
# purple_cnv_files <- (function(){
#    path <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/ploidy_estimator/purple_cnv_paths.rds'
#    
#    if(!file.exists(path)){
#       purple_cnv_files <- list.files(pattern="\\.purple.cnv$", path="/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-DR047/data/", recursive=T, full.names=T)
#       saveRDS(purple_cnv_files, path)
#    } else {
#       readRDS(path)
#    }
# })()
# 
# 
# pb <- txtProgressBar(max=length(purple_cnv_files), style=3)
# counter <- 0
# arm_ploidies <- lapply(purple_cnv_files, function(i){
#    counter <<- counter + 1
#    setTxtProgressBar(pb, counter)
#    calcChromArmPloidies(i)
# })
# 
# # chrom_arm_names <- paste0(rep(c(1:22,'X'),each=2),c('p','q'))
# # chrom_arm_names <- chrom_arm_names[!(chrom_arm_names %in% paste0(c(13,14,15,21,22),'p'))] ## Keep only q arm for one arm chromosomes
# # chrom_arm_names <- c(chrom_arm_names,'genome')
# m2 <- do.call(rbind,lapply(arm_ploidies, function(i){
#    i[chrom_arm_names]
# }))
# rownames(m2) <- sapply(strsplit(basename(purple_cnv_files),'[.]'), `[[`, 1)
# m2 <- as.data.frame(m2)
# write.table(
#    m2, '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/ploidy_estimator/arm_ploidies.txt',
#    quote=F, sep='\t'
# )
# 
# ###
# m1 <- read.delim('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/ploidy_estimator/melanie/Results_CopyNumber_withX.tsv')
# m1 <- t(m1)
# colnames(m1) <- m1[1,]; m1 <- m1[-1,]
# sample_names <- rownames(m1)
# m1 <- apply(m1,2,as.numeric)
# rownames(m1) <- sample_names
# m1 <- as.data.frame(m1)
# m1 <- m1[!apply(m1,1,anyNA),]
# ###
# 
# library(Rtsne)
# metadata <- read.delim('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_DR010_DR047/metadata/DR-047-update1_metadata.txt')
# 
# getMetadata <- function(v,key='primary_tumor_location'){
#    #v=rownames(sbs_sigs)
#    metadata_match <- match(v,metadata$sample_id)
#    metadata[metadata_match,key]
# }
# 
# mkTsnePd <- function(m, seed=1, ...){
#    #m <- sbs_sigs/rowSums(sbs_sigs)
#    set.seed(seed)
#    tsne <- Rtsne(m, dims=2, check_duplicates=F, verbose=T, ...)
#    
#    df <- as.data.frame(tsne$Y)
#    colnames(df) <- c('x','y')
#    df$sample <- rownames(m) 
#    df$cancer_type <- getMetadata(df$sample)
#    
#    return(df)
# }
# 
# CANCER_TYPE_COLORS <- c(
#    "Breast"='#E41A1C',
#    "Colon/Rectum"='#377EB8',
#    "Lung"='#4DAF4A',
#    "Prostate"='#984EA3',
#    "Skin"='#FF7F00',
#    "Bone/Soft tissue"='#FFFF33',
#    "Esophagus"='#A65628',
#    "Urinary tract"='#F781BF',
#    "Ovary"='#999999',
#    "Kidney"='#8DD3C7',
#    "NET"='#FFFCBB',
#    "Biliary"='#BEBADA',
#    
#    "Unknown"='#F0F0F0',
#    
#    "Pancreas"='#FB8072',
#    "Nervous system"='#80B1D3',
#    "Uterus"='#FDB462',
#    "Head and neck"='#B3DE69',
#    "Liver"='#FCCDE5',
#    "Stomach"='#D9D9D9',
#    "Mesothelioma"='#BC80BD',
#    "Lymphoid"='#CCEBC5',
#    "Thyroid"='#912643',
#    "Small intestine"='#143167',
#    "Vulva"='#194320',
#    "Penile"='#461348',
#    "Testis"='#4F310F',
#    "Adrenal"='#B5367B',
#    "Thymus"='#B4D8E7',
#    "Myeloid"='#F5CFB0',
#    "Double primary"='#D54088',
#    "Eye"='black'
# )
# 
# plotTsne <- function(pd, plot.title=NULL, out.path=NULL){
#    p <- ggplot(pd,aes(x,y,color=cancer_type)) +
#       geom_point() +
#       labs(x='t-SNE dim. 2', y='t-SNE dim. 1', color='Cancer type') +
#       scale_color_manual(values=CANCER_TYPE_COLORS) +
#       theme_bw() +
#       theme(
#          plot.title=element_text(hjust=0.5)
#       )
#    
#    if(!is.null(plot.title)){ p <- p + ggtitle(plot.title) }
#    
#    if(!is.null(out.path)){ 
#       pdf(out.path,10,7)
#       grid.draw(p)
#       dev.off()
#    } else {
#       return(p)
#    }
# }
# 
# pd1 <- mkTsnePd(m1)
# plotTsne(pd1)
# 
# pd2 <- mkTsnePd(m2)
# plotTsne(pd2)
# 
# table(m2$genome)
# 
# m2[m2$genome==1,]
# 
# m2[order(m2$genome, decreasing=T),]



