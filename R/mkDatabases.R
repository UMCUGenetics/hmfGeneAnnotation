####################################################################################################
#' Download ClinVar vcf and convert to txt file 
#'
#' @param out.file Path to output file
#' @param tmp.dir Temporary processing directory. Defaults to ~/
#' @param java.path Path to java binary (defaults to the installed JRE location)
#' @param snpsift.path Path to SnpSift jar (defaults to the one included in this package)
#' @param verbose Show progress messages?
#'
#' @export
#'
#' @example
#' mkClinvarTxt('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/db/clinvar.txt.bgz')
mkClinvarTxt <- function(
   out.file, tmp.dir='~/',
   url='ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20200224.vcf.gz',
   java.path=JAVA_PATH, snpsift.path=SNPSIFT_PATH,
   verbose=T
){
   #out.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/db/clinvar.txt.bgz'
   
   tmp_paths <- list()
   tmp_paths$clinvar_vcf <- paste0(tmp.dir,'/clinvar.vcf.gz')
   
   if(verbose){ message('Downloading clinvar vcf to ',tmp.dir) }
   if(!file.exists(tmp_paths$clinvar_vcf)){
      download.file(url,tmp_paths$clinvar_vcf)
   }
   
   if(verbose){ message('Extracting relevant data from vcf to txt...') }
   tmp_paths$clinvar_txt <- paste0(tmp.dir,'/clinvar.txt.gz')
   fileConn <- gzfile(tmp_paths$clinvar_txt)
   writeLines(
      '#chrom\tpos\tref\talt\tsig\tid', 
      fileConn
   )
   close(fileConn)
   
   string <- paste0(
      'JAVA=',java.path,'\n',
      'SNPSIFT=',snpsift.path,'\n',
      'VCF_FILE=',tmp_paths$clinvar_vcf,'\n',
      'OUT_FILE=',tmp_paths$clinvar_txt,'\n',
      
      '$JAVA -jar $SNPSIFT extractFields -v $VCF_FILE ',
      'CHROM POS REF ALT CLNSIG[0] ID | \n',
      'tail -n+2 | gzip -c >> $OUT_FILE'
   )
   
   system(string)
   #file.remove(clinvar_vcf_tmp_path)
   
   if(verbose){ message('bgzipping and making tabix index...') }
   
   if(file.exists(out.file)){
      warning('Removed existing output file:', out.file)
      file.remove(out.file)
   }
   
   Rsamtools::bgzip(tmp_paths$clinvar_txt, dest=out.file)
   Rsamtools::indexTabix(out.file, seq=1, start=2, end=2, skip=1)
   
   if(verbose){ message('Removing tmp files...') }
   for(i in tmp_paths){ file.remove(i) }
}

#clinvar_txt <- read.delim('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/db/clinvar.txt.bgz', stringsAsFactors=F)


####################################################################################################
#' Retrieve gene identifiers from HGNC
#'
#' @param hgnc.url URL to HGNC database, where the output contains the columns: HGNC ID, Approved
#' symbol, Previous symbols, Synonyms, ENSEMBL gene id
#' @param out.path Export path
#' @param export.as.rdata Save as RData file?
#' @param verbose Show messages?
#' @param ... Arguments that can be passed to biomaRt::useMart()
#'
#' @return A dataframe of gene identifers from HGNC
#' @export
#'
retrieveHgncGeneList <- function(
   ## excludes withdrawn symbols
   hgnc.url='https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit',
   out.path=NULL,
   export.as.rdata=F,
   verbose=T,
   ...
){
   if(verbose){ message('Downloading HGNC gene list...') }
   genes_hgnc <- read.delim(hgnc.url,check.names=F, comment.char='#',stringsAsFactors=F)
   colnames(genes_hgnc) <- gsub(' ','_',tolower(colnames(genes_hgnc)))
   
   mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org',...)
   #mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org')
   
   empty_ensg <- nchar(genes_hgnc$ensembl_gene_id) == 0
   if(any(empty_ensg) & verbose){
      message(sprintf(
         '%s/%s entries were found without ENSGs. Attempting to retrieve ENSGs using biomaRt...', 
         sum(empty_ensg), nrow(genes_hgnc)
      ))
   }
   
   ## Split dataframe into missing/non-missing ENSGs
   genes_hgnc_with_ensg <- genes_hgnc[!empty_ensg,]
   genes_hgnc_no_ensg <- genes_hgnc[empty_ensg,]
   
   ## Retrieve ENSGs
   biomart_out <- getBM(
      mart=useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org'),
      attributes=c('hgnc_symbol','ensembl_gene_id'),
      filters='hgnc_symbol',
      values=genes_hgnc_no_ensg$approved_symbol,
      verbose=F
   )
   
   new_ensg <- biomart_out$ensembl_gene_id[match(genes_hgnc_no_ensg$approved_symbol, biomart_out$hgnc_symbol)]
   
   if(any(is.na(new_ensg)) & verbose){
      message(sprintf(
         'ENSGs for %s/%s entries could be retrieved using biomaRt...', 
         sum(!is.na(new_ensg)), length(new_ensg)
      ))
   }
   
   ## Return
   new_ensg[is.na(new_ensg)] <- ''
   genes_hgnc_no_ensg$ensembl_gene_id <- new_ensg
   
   GENES_HGNC <- rbind(genes_hgnc_with_ensg, genes_hgnc_no_ensg)
   GENES_HGNC <- GENES_HGNC[match(GENES_HGNC$hgnc_id, genes_hgnc$hgnc_id),]
   
   if(is.null(out.path)){
      return(GENES_HGNC)
   } else {
      message(sprintf('Exporting table to: %s', out.path))
      if(export.as.rdata & verbose){
         save(GENES_HGNC, file=out.path)
      } else {
         write.table(GENES_HGNC, gzfile(out.path), sep='\t',row.names=F,quote=F)
      }
   }
}


####################################################################################################
mkGenesBed <- function(
   cosmic.genes.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/db/cosmic_cancer_gene_census_20200225.tsv.gz',
   out.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/misc/cosmic_cancer_gene_census_20200225.bed'
){

   cosmic_genes <- read.delim(cosmic.genes.path, stringsAsFactors=F)
   colnames(cosmic_genes) <- gsub('[.]','_',colnames(cosmic_genes))
   colnames(cosmic_genes) <- gsub('_$','',colnames(cosmic_genes))
   colnames(cosmic_genes) <- tolower(colnames(cosmic_genes))
   
   ensg_ids <- regmatches(cosmic_genes$synonyms, gregexpr('ENSG\\d+[.]', cosmic_genes$synonyms))
   ensg_ids <- unlist(lapply(ensg_ids, function(i){
      if(length(i)==0){ NA } else { i }
   }))
   ensg_ids <- gsub('[.]','',ensg_ids)
   
   cosmic_genes$ensembl_gene_id <- ensg_ids
   cosmic_genes$hgnc_symbol <- ensgToHgncSymbol(ensg_ids)
   cosmic_genes_ss <- cosmic_genes[!is.na(cosmic_genes$hgnc_symbol),]
   
   mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org')
   
   genes_bed <- getBM(
      mart=mart, verbose=F,
      attributes=c(
         'chromosome_name', 'start_position', 'end_position',
         'hgnc_id','hgnc_symbol','ensembl_gene_id'
      ),
      filters='ensembl_gene_id',
      values=ensg_ids
   )
   
   ## Clean up table
   common_ensg_ids <- intersect(cosmic_genes_ss$ensembl_gene_id, genes_bed$ensembl_gene_id)
   genes_bed <- genes_bed[
      genes_bed$ensembl_gene_id %in% common_ensg_ids & ## Only keep genes with ENSG
      !apply(genes_bed,1,anyNA) & ## Rm rows with NA
      !grepl('PATCH',genes_bed$chromosome_name) ## Rm genome patch genes
   ,]
   
   ## Verify that ENSG and HGNC ids can be converted back and forth
   # anyNA(ensgToHgncSymbol(genes_bed$ensembl_gene_id))
   # anyNA(geneNamesToEnsg(genes_bed$hgnc_symbol))
   
   ## Assign TSG or oncogene
   genes_bed$role_in_cancer <- cosmic_genes$role_in_cancer[ match(genes_bed$ensembl_gene_id, cosmic_genes$ensembl_gene_id) ]
   genes_bed$role_in_cancer <- gsub(' ','',genes_bed$role_in_cancer)
   genes_bed$role_in_cancer[nchar(genes_bed$role_in_cancer)==0] <- 'unknown'
   
   genes_bed$molecular_genetics <- cosmic_genes$molecular_genetics[ match(genes_bed$ensembl_gene_id, cosmic_genes$ensembl_gene_id) ]
   genes_bed$molecular_genetics[nchar(genes_bed$molecular_genetics)==0] <- 'unknown'
   
   genes_bed <- genes_bed[order(genes_bed$hgnc_symbol),]
   
   colnames(genes_bed)[1:3] <- c('#chrom','start','end')
   
   write.table(genes_bed, out.path, sep='\t', row.names=F, quote=F)
}

####################################################################################################
getCentromerePositions <- function(
   out.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/db/centromere_positions_hg19.txt'
){
   ## Downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/gap.txt.gz
   gap <- read.table(
      '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/db/gap.txt.gz',
      stringsAsFactors=F
   )
   gap <- gap[,c(2:4,8)]
   colnames(gap) <- c('chrom','start','end','gap_type')
   
   out <- subset(gap, gap_type=='centromere')
   out$pos <- (out$start + out$end)/2
   out <- out[naturalsort::naturalorder(out$chrom),c('chrom','pos')]
   out$chrom <- gsub('chr','',out$chrom)
   
   write.table(out, out.path, sep='\t', row.names=F, quote=F)
}

# getChromArmCoords <- function(
#    out.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/db/chrom_arms_hg19.txt'
# ){
#    ## Not how HMF calculates chrom arms
#    ## Downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz
#    cytobands <- read.table(
#       '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/db/cytoBand.txt.gz',
#       stringsAsFactors=F
#    )
#    colnames(cytobands) <- c('chrom','start','end','cytoband','stain_value')
#    cytobands$arm <- substring(cytobands$cytoband,1,1)
#    cytobands_by_chrom <- split(cytobands, cytobands$chrom)
# 
#    arm_coords_by_chrom <- lapply(cytobands_by_chrom, function(i){
#       #i=cytobands_by_chrom$chr21
#       df <- do.call(rbind,lapply(split(i, i$arm), function(i){
#          ## Select 1st and last rows
#          data.frame(
#             chrom=i$chrom[1],
#             arm=i$arm[1],
#             start=i$start[1]+1, ## Make coords 1 based
#             end=i$end[nrow(i)]
#          )
#       }))
# 
#       if(df$chrom[1] %in% paste0('chr',ONE_ARMED_CHROMS)){
#          df <- data.frame(
#             chrom=df$chrom[1],
#             arm='q',
#             start=i$start[1],
#             end=i$end[nrow(i)]
#          )
#       }
# 
#       return(df)
#    })
# 
#    arm_coords_by_chrom <- arm_coords_by_chrom[naturalsort::naturalorder(names(arm_coords_by_chrom))]
#    out <- do.call(rbind, arm_coords_by_chrom)
#    rownames(out) <- NULL
# 
#    write.table(out, out.path, sep='\t', row.names=F, quote=F)
# }










