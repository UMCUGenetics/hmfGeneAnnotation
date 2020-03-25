####################################################################################################
#' Subsets PURPLE gene cnv file for genes of interest
#'
#' @param gene.cnv.file Path to PURPLE gene cnv file
#' @param out.file Path to output file
#' @param bed.file Path to bed file specifying the genome intervals to keep (defaults ot the one 
#' stored in this package)
#' @param genes.enst2ensg.file Path to lookup table to convert ENSEMBL transcript to genes ids
#' @param verbose Show progress messages?
#'
#' @return Nothing but writes the subsetted PURPLE gene cnv txt file
#' @export
#'
preProcessGeneCnv <- function(
   gene.cnv.file, out.file, 
   bed.file=BED_FILE, 
   genes.enst2ensg.file=GENES_ENST2ENSG,
   verbose=T
){
   
   # gene.cnv.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR-104/data/somatics/191205_HMFregCPCT_FR16670564_FR17496616_CPCT02020989/CPCT02020989T.purple.cnv.gene.tsv'
   # out.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02010349T.purple.gene.cnv.txt'
   
   if(verbose){ message('Reading lookup tables...') }
   genes_enst2ensg <- read.delim(genes.enst2ensg.file, stringsAsFactors=F, check.names=F)
   genes_bed <- read.delim(bed.file, stringsAsFactors=F, check.names=F)
   
   if(verbose){ message('Reading gene cnv file...') }
   gene_cnv <- read.delim(gene.cnv.file, stringsAsFactors=F)
   
   ## Convert colnames from camelcase to snakecase
   colnames(gene_cnv) <- gsub("([a-z])([A-Z])", "\\1_\\L\\2", colnames(gene_cnv), perl = TRUE)
   colnames(gene_cnv) <- tolower(colnames(gene_cnv))
   
   ## Retrieve ensembl gene ids
   if(verbose){ message('Retrieving ensembl gene ids...') }

   ins_cols <- genes_enst2ensg[
      match(gene_cnv$transcript_id, genes_enst2ensg$ensembl_transcript_id),
      c('ensembl_gene_id','hgnc_symbol')
   ]
   gene_cnv <- cbind(gene_cnv,ins_cols)
   
   ## Subset by ensembl gene id of selected genes
   if(verbose){ message('Subsetting gene cnv table...') }
   gene_cnv <- gene_cnv[gene_cnv$ensembl_gene_id %in% genes_bed$ensembl_gene_id,]
   
   ## Remove duplicate gene ids
   gene_cnv <- gene_cnv[!(duplicated(gene_cnv$ensembl_gene_id)),]
   
   ## Determine chrom arm
   gene_cnv$chrom_arm <- paste0(
      gene_cnv$chromosome,
      substring(gene_cnv$chromosome_band,1,1)
   )
   
   ## Select relevant rows
   sel_cols <- c(
      'chromosome','start','end','chrom_arm','hgnc_symbol','ensembl_gene_id',
      'min_copy_number','max_copy_number','min_minor_allele_ploidy'
   )
   gene_cnv <- gene_cnv[,sel_cols]
   colnames(gene_cnv)[1] <- 'chrom'
   
   ## Export
   if(verbose){ message('Exporting new gene cnv table...') }
   write.table(gene_cnv, out.file, sep='\t', row.names=F, quote=F)
}

# subsetGeneCnv(
#    gene.cnv.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR-104/data/somatics/191205_HMFregCPCT_FR16670564_FR17496616_CPCT02020989/CPCT02020989T.purple.cnv.gene.tsv',
#    out.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02010349T.purple.gene.cnv.txt'
# )

####################################################################################################
#' Filter vcf for compatibility with downstream functions
#' 
#' @description Wrapper function for SnpSift filter and SnpSift intervals. Selects PASS variants and 
#' variants at coords specified in bed file. For germline vcfs, remove somatic variants.
#' 
#' @param vcf.file Path to vcf file (gzip compressed)
#' @param out.file Path to output vcf file (make sure to add .gz at the end)
#' @param mode Can be 'germ' or 'som'
#' @param bed.file Path to bed file specifying the genome intervals to keep (defaults ot the one 
#' stored in this package)
#' @param java.path Path to java binary (defaults to the installed JRE location)
#' @param snpsift.path Path to SnpSift jar (defaults to the one included in this package)
#'
#' @return Nothing but writes a gzipped vcf file
#' @export
#'
filterVcf <- function(
   vcf.file, out.file, mode=NULL, bed.file=BED_FILE,
   java.path=JAVA_PATH, snpsift.path=SNPSIFT_PATH
){
   #vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR-104/data/somatics/191205_HMFregCPCT_FR16670564_FR17496616_CPCT02020989/CPCT02020989T.purple.somatic.vcf.gz'
   #out.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02020989T.purple.somatic.ss.vcf.gz'
   
   if(mode=='som'){
      filter_string <- "\"(FILTER='PASS')\""
   } else if(mode=='germ'){
      ## Remove somatic variants (i.e GT with 0/0 or ./. in the germline (1st) column)
      filter_string <- "\"(FILTER='PASS') & !(GEN[0].GT='0/0') & !(GEN[0].GT='./.')\"" 
   } else {
      message("Mode unspecified. Filtering intervals only")
   }
   
   string <- paste0(
      'JAVA=',java.path,'\n',
      'SNPSIFT=',snpsift.path,'\n',
      'VCF_FILE=',vcf.file,'\n',
      'BED_FILE=',bed.file,'\n',
      'OUT_FILE=',out.file,'\n',
      'FILTER_STRING=',filter_string,'\n',
      
      'gunzip -c $VCF_FILE | \n',
      if(!is.null(mode)){ '$JAVA -jar $SNPSIFT filter -v "$FILTER_STRING" | \n' } else { '' },
      '$JAVA -jar $SNPSIFT intervals -verbose $BED_FILE | \n',
      'gzip -c > $OUT_FILE'
   )
   
   # fileConn<-file("/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/test.sh")
   # writeLines(string, fileConn)
   # close(fileConn)
   
   system(string)
}
# 
# filterVcf(
#    vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR-104/data/somatics/191205_HMFregCPCT_FR16670564_FR17496616_CPCT02020989/CPCT02020989T.purple.somatic.vcf.gz',
#    out.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02020989T.purple.somatic.ss.vcf.gz',
#    mode='som'
# )

#---------------------------------------------------------------------------------------------------
#' Annotate variant type with snpEff
#'
#' @param vcf.file Path to vcf file (gzip compressed)
#' @param out.file Path to output txt file (make sure to add .gz at the end)
#' @param genome Name of the genome to supply to snpEff
#' @param java.path Path to java binary (defaults to the installed JRE location)
#' @param snpsift.path Path to SnpSift jar (defaults to the one included in this package)
#'
#' @return Nothing but writes a gzipped txt file
#' @export
#'
annotateVariantType <- function(
   vcf.file, out.file, genome='GRCh37.75',
   java.path=JAVA_PATH, snpeff.path=SNPEFF_PATH
){
   #vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR-104/data/somatics/191205_HMFregCPCT_FR16670564_FR17496616_CPCT02020989/CPCT02020989T.purple.somatic.vcf.gz'
   #out.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02020989T.purple.somatic.ss.vcf.gz'
   
   string <- paste0(
      'JAVA=',java.path,'\n',
      'SNPEFF=',snpeff.path,'\n',
      'GENOME=',genome,'\n',
      'VCF_FILE=',vcf.file,'\n',
      'OUT_FILE=',out.file,'\n',

      '$JAVA -jar $SNPEFF ann $GENOME $VCF_FILE -lof -no-downstream -no-intergenic -noShiftHgvs -verbose | \n',
      'gzip -c > ${OUT_FILE}.tmp \n',
      'mv ${OUT_FILE}.tmp ${OUT_FILE}' ## Hack to do annotation inline of vcf.file is same as out.file
   )
   
   system(string)
}

# annotateVariantType(
#    vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02020989T.purple.somatic.ss.vcf.gz',
#    out.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02020989T.purple.somatic.ss2.vcf.gz'
# )

#---------------------------------------------------------------------------------------------------
#' Extracts revelant fields for downstream functions and outputs as txt file
#'
#' @param vcf.file Path to vcf file (gzip compressed)
#' @param out.file Path to output txt file (make sure to add .gz at the end)
#' @param java.path Path to java binary (defaults to the installed JRE location)
#' @param snpsift.path Path to SnpSift jar (defaults to the one included in this package)
#'
#' @return Nothing but writes a gzipped txt file
#' @export
#'
extractVcfFields <- function(
   vcf.file, out.file,
   java.path=JAVA_PATH, snpsift.path=SNPSIFT_PATH
){
   # vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02020989T.purple.somatic.ss.vcf.gz'
   # out.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02020989T.purple.somatic.ss.txt.gz'
   
   fileConn <- gzfile(out.file)
   writeLines(
      'chrom\tpos\tref\talt\tsnpeff_eff\tsnpeff_gene\tensembl_gene_id\thgvs_c', 
      fileConn
   )
   close(fileConn)
   
   string <- paste0(
      'JAVA=',java.path,'\n',
      'SNPSIFT=',snpsift.path,'\n',
      'VCF_FILE=',vcf.file,'\n',
      'OUT_FILE=',out.file,'\n',
      
      '$JAVA -jar $SNPSIFT extractFields -v $VCF_FILE ',
      'CHROM POS REF ALT ANN[0].EFFECT ANN[0].GENE ANN[0].GENEID ANN[0].HGVS_C | \n',
      'tail -n+2 | gzip -c >> $OUT_FILE'
   )
   
   system(string)
}

# extractVcfFields(
#    vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02020989T.purple.somatic.ss.vcf.gz',
#    out.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts_prototype/test_output/CPCT02020989T.purple.somatic.ss.txt.gz'
# )


####################################################################################################
#' Get clinical significance by querying a tabix indexed bed file
#'
#' @param df A bed file like dataframe containing the columns: chrom, pos, ref, alt
#' @param db.path Path to a tabix indexed bed file with the clinical significance annotations (e.g.
#' from ClinVar). This bed file should contain the columns (in this order): chrom, pos, ref, alt, 
#' clinsig, id. However, it should not have colnames.
#' @param verbose Show progress messages?
#'
#' @return A character vector of the clinical significance annotations
#' @export
#'
getClinSig <- function(df, db.path, verbose=T){

   ##
   if(verbose){ message('Querying database using coords from df...') }
   coords <- paste0(df$chrom,':',df$pos,'-',df$pos)
   #coords <- '1:1-1'
   
   tabix_out <- seqminer::tabix.read(db.path, coords)
   
   if(length(tabix_out)==0){ return(rep(NA,nrow(df))) }
   
   tabix_out <- lapply(strsplit(tabix_out,'\t'), as.vector)
   tabix_out <- as.data.frame(do.call(rbind,tabix_out), stringsAsFactors=F)
   colnames(tabix_out) <- c('chrom','pos','ref','alt','clinsig','id')
   tabix_out$pos <- as.integer(tabix_out$pos)
   
   ##
   if(verbose){ message('Separating df rows with/without entry in database...') }
   tabix_coords <- paste0(tabix_out$chrom,':',tabix_out$pos,'-',tabix_out$pos)
   df_split <- split(df, coords %in% tabix_coords)
   names(df_split) <- c(
      'na', ## These rows don't have an entry in the tabix output
      'proc' ## These rows have the coords, but need to check for ref and alt match
   )
   
   ##
   if(verbose){ message('Retrieving clinical significance for df rows with entry in database...') }
   clinsigs <- with(df_split$proc, {
      unlist(Map(function(chrom,pos,ref,alt){
         
         query <- tabix_out[
            tabix_out$chrom==chrom & tabix_out$pos==pos & tabix_out$ref==ref & tabix_out$alt==alt,
            'clinsig'
            ]
         
         if(length(query)==0){ NA } else { query }
         
      },chrom,pos,ref,alt))
   })
   
   clinsigs <- c(
      rep(NA,nrow(df_split$na)), ## For rows without entry in database, return NA
      clinsigs
   )
   
   ## Restore order of original df
   if(verbose){ message('Returning output...') }
   row_order <- as.integer(c(rownames(df_split$na), rownames(df_split$proc)))
   clinsigs[order(row_order)]

}

# df <- read.delim(
#   '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CUPs_classifier_data/gene_ann/CPCT02010349T/proc/CPCT02010349T.germ.txt.gz',
#   stringsAsFactors=F
# )
# 
# df$clinsig <- getClinSig(
#   df,
#   db.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/db/clinvar.txt.bgz'
# )
# 
# View(df[!is.na(df$clinsig),])
