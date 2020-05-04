####################################################################################################
#' Determine chromosome arm
#'
#' @param chrom A character vector specifying the chromosomes
#' @param pos An integer vector specifying the positions
#' @param centromere.positions.path Path to table containing for each chrom the centromere positions,
#' and whether the chrom is one armed. Table should have the columns chrom, pos (centromere pos),
#' one_armed 
#'
#' @return
#' @export
#'
#' @examples
getChromArms <- function(chrom, pos, centromere.positions.path=CENTROMERE_POSITIONS){
   #chrom=gene_cnv$chrom
   #pos=gene_cnv$start
   #v=CENTROMERE_POSITIONS
   
   centromere_positions <- read.delim(centromere.positions.path, stringsAsFactors=F)
   
   if(length(chrom)!=length(pos)){
      stop('chrom and pos are not the same length')
   }
   
   df <- data.frame(chrom,pos)
   
   df$centro_pos <- centromere_positions[match(df$chrom, centromere_positions$chrom),'pos']
   df$one_armed <- centromere_positions[match(df$chrom, centromere_positions$chrom),'one_armed']
   
   df$arm <- 'p'
   df$arm[df$pos>=df$centro_pos] <- 'q'
   df$arm[df$one_armed] <- 'q'
   
   paste0(df$chrom, df$arm)
}


#---------------------------------------------------------------------------------------------------
#' Extract gene copy number info from cnv table
#'
#' @param cnv.file A copy number txt file with the columns chrom, start, end, total_cn, major_cn (in
#' that order)
#' @param exons.bed.file Path to the bed file containing the columns chrom, exon_start, exon_end, 
#' hgnc_symbol, ensembl_gene_id
#' @param out.file If not NULL will write the output txt to this path
#'
#' @details Was originally written to deal with PCAWG CNV files
#'
#' @return Writes a txt file or returns a dataframe
#' @export
#'
cnvToGeneCnv <- function(cnv.file, exons.bed.file=EXONS_BED_FILE, out.file=NULL, verbose=T){
   #cnv.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/PCAWG_2020/vcf/somatic/cnv/cna_annotated/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20170119.somatic.cna.annotated.txt'
   #exons.bed.file=EXONS_BED_FILE

   if(verbose){ message('Reading cnv file...') }
   cnv <- read.delim(cnv.file, stringsAsFactors=F)
   col_names <- c('chrom','start','end','total_cn','major_cn','minor_cn')
   colnames(cnv) <- col_names
   cnv <- cnv[,col_names]

   if(verbose){ message('Reading bed file...') }
   bed <- read.delim(exons.bed.file, stringsAsFactors=F)
   colnames(bed)[1] <- 'chrom'

   if(verbose){ message('Finding CNVs that overlap with bed file exons...') }
   exon_cnv <- lapply(1:nrow(bed), function(i){
      #i=9937
      row <- bed[i,]

      ## Get any cnv interval that overlaps with exon interval
      df <- cnv[

         cnv$chrom==row$chrom

         & (
            ## start or end of cnv is between gene
            (row$exon_start<=cnv$start & cnv$start<=row$exon_end)
            | (row$exon_start<=cnv$end & cnv$end<=row$exon_end)

            ## exon is entirely within cnv
            | (cnv$start<=row$exon_start & cnv$end>=row$exon_end)
         )
      ,]

      out <- row[,c('chrom','exon_start','exon_end','hgnc_symbol','ensembl_gene_id')]
      rownames(out) <- NULL

      if(nrow(df)==0){
         values <- data.frame(
            total_cn=NA,
            minor_cn=NA,
            row.names=NULL
         )
      } else {
         values <- data.frame(
            total_cn=df$total_cn,
            minor_cn=df$minor_cn,
            row.names=NULL
         )
      }

      # if(nrow(df)==0){
      #    values <- data.frame(
      #       min_copy_number=NA,
      #       max_copy_number=NA,
      #       min_minor_allele_ploidy=NA
      #    )
      # } else if(nrow(df)==1){
      #    values <- data.frame(
      #       min_copy_number=df$total_cn,
      #       max_copy_number=df$total_cn,
      #       min_minor_allele_ploidy=df$minor_cn
      #    )
      # } else {
      #    values <- data.frame(
      #       min_copy_number=min(df$total_cn, na.rm=T),
      #       max_copy_number=max(df$total_cn, na.rm=T),
      #       min_minor_allele_ploidy=min(df$minor_cn, na.rm=T)
      #    )
      # }

      cbind(out, values)
   })
   exon_cnv <- do.call(rbind, exon_cnv)

   if(verbose){ message('Calculating CN info per gene...') }
   exon_cnv_split <- split(exon_cnv, exon_cnv$ensembl_gene_id)
   gene_cnv <- do.call(rbind, lapply(exon_cnv_split, function(i){
      data.frame(
         ensembl_gene_id=i[1,'ensembl_gene_id'],
         min_copy_number=min(i$total_cn, na.rm=T),
         max_copy_number=max(i$total_cn, na.rm=T),
         min_minor_allele_ploidy=min(i$minor_cn, na.rm=T)
      )
   }))
   rownames(gene_cnv) <- NULL
   
   gene_cnv$min_copy_number[!is.finite(gene_cnv$min_copy_number)] <- NA
   gene_cnv$max_copy_number[!is.finite(gene_cnv$max_copy_number)] <- NA
   gene_cnv$min_minor_allele_ploidy[!is.finite(gene_cnv$min_minor_allele_ploidy)] <- NA

   gene_info <- unique(bed[,c('chrom','start','end','hgnc_symbol','ensembl_gene_id')])
   gene_cnv <- cbind(
      gene_info,
      gene_cnv[match(gene_info$ensembl_gene_id, gene_cnv$ensembl_gene_id),-1]
   )

   if(verbose){ message('Annotating chrom arm...') }
   gene_cnv$chrom_arm <- getChromArms(gene_cnv$chrom, gene_cnv$start)

   if(!is.null(out.file)){
      if(verbose){ message('Writing output...') }
      write.table(gene_cnv, out.file, sep='\t', row.names=F, quote=F)
   } else {
      return(gene_cnv)
   }
}

# cnvToGeneCnv <- function(cnv.file, bed.file=BED_FILE, out.file=NULL){
#    #cnv.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/PCAWG_2020/vcf/somatic/cnv/cna_annotated/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20170119.somatic.cna.annotated.txt'
#    #cnv.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/PCAWG_2020/vcf/cnv/cna_annotated/43b675e0-22e0-42d6-a060-afc93e22ac95.consensus.20170119.somatic.cna.annotated.txt'
#    #bed.file=BED_FILE
#    
#    if(verbose){ message('Reading cnv file...') }
#    cnv <- read.delim(cnv.file)
#    col_names <- c('chrom','start','end','total_cn','major_cn','minor_cn')
#    colnames(cnv) <- col_names
#    cnv <- cnv[,col_names]
#    
#    if(verbose){ message('Reading bed file...') }
#    bed <- read.delim(bed.file)
#    colnames(bed)[1] <- 'chrom'
#    
#    if(verbose){ message('Finding CNVs that overlap with bed file intervals...') }
#    gene_cnv <- lapply(1:nrow(bed), function(i){
#       #i=230
#       gene <- bed[i,]
#       
#       ## Get any cnv interval that overlaps with gene interval
#       df <- cnv[
#          
#          cnv$chrom==gene$chrom
#          
#          & (
#             ## start or end of cnv is between gene
#             (gene$start<=cnv$start & cnv$start<=gene$end)
#             | (gene$start<=cnv$end & cnv$end<=gene$end)
#             
#             ## gene is entirely within cnv
#             | (cnv$start<=gene$start & cnv$end>=gene$end)
#          )
#          ,]
#       
#       out <- gene[,c('chrom','start','end','hgnc_symbol','ensembl_gene_id')]
#       
#       if(nrow(df)==0){
#          values <- data.frame(
#             min_copy_number=NA,
#             max_copy_number=NA,
#             min_minor_allele_ploidy=NA
#          )
#       } else if(nrow(df)==1){
#          values <- data.frame(
#             min_copy_number=df$total_cn,
#             max_copy_number=df$total_cn,
#             min_minor_allele_ploidy=df$minor_cn
#          )
#       } else {
#          values <- data.frame(
#             min_copy_number=min(df$total_cn, na.rm=T),
#             max_copy_number=max(df$total_cn, na.rm=T),
#             min_minor_allele_ploidy=min(df$minor_cn, na.rm=T)
#          )
#       }
#       
#       cbind(out, values)
#    })
#    gene_cnv <- do.call(rbind, gene_cnv)
#    
#    if(verbose){ message('Annotating chrom arm...') }
#    gene_cnv$chrom_arm <- getChromArms(gene_cnv$chrom, gene_cnv$start)
#    
#    if(!is.null(out.file)){
#       if(verbose){ message('Writing output...') }
#       write.table(gene_cnv, out.file, sep='\t', row.names=F, quote=F)
#    } else {
#       return(gene_cnv)
#    }
# }

#---------------------------------------------------------------------------------------------------
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
  
   ## Select relevant rows
   sel_cols <- c(
      'chromosome','start','end','hgnc_symbol','ensembl_gene_id',
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
      filter_string <- "\"(na FILTER) | (FILTER='PASS')\""
   } else if(mode=='germ'){
      ## Remove somatic variants (i.e GT with 0/0 or ./. in the germline (1st) column)
      filter_string <- "\"((na FILTER) | (FILTER='PASS')) & !(GEN[0].GT='0/0') & !(GEN[0].GT='./.')\"" 
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

      '$JAVA -jar $SNPEFF ann $GENOME $VCF_FILE -lof -no-downstream -no-intergenic -noShiftHgvs -noStats -verbose | \n',
      'gzip -c > ${OUT_FILE}.tmp \n',
      'mv ${OUT_FILE}.tmp ${OUT_FILE}' ## Hack to do annotation inline of vcf.file
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
