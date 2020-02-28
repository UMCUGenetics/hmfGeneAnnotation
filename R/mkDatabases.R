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
   
   genes_bed <- genes_bed[
      naturalorder(paste0(genes_bed$chromosome_name,'_',genes_bed$start_position))
   ,]
   
   colnames(genes_bed)[1:3] <- c('#chrom','start','end')
   
   write.table(genes_bed, out.path, sep='\t', row.names=F, quote=F)
}
