library(magrittr)
library(stringr)
library(reshape2)
library(ggplot2)
library(scales)
library(RColorBrewer)
options(stringsAsFactors = F)

#========= Paths =========#
base_dir <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'
if(dir.exists('/Users/lnguyen/')){
   base_dir <- paste0('/Users/lnguyen/', base_dir)
}

#========= Load data =========#
in_dir <- paste0(base_dir,'/scripts_main/hmfVariantAnnotation/analysis/hr_gene_def/')
gene_statuses <- read.table(paste0(base_dir,'/HMF_update/analysis/hr_gene_def/hmf_gene_statuses.txt'),sep='\t',header=T)

## CHORD output
chord_pred_path <- paste0(base_dir,'/scripts_main/mltoolkit/random_forest_training_v2/models/37_HMF_pkgMutSigs_germLt3-somEq0-lohEq0_greylistIncl_multiclass_snvContext_svContext_allSigsRel_scaleGridSearch_featPreFilt_boruta_newAnn/whole_hmf_dataset_probs.txt')
chord_pred <- read.table(chord_pred_path,sep='\t',header=T)

## metadata
metadata_path <- paste0(base_dir,'/HMF_update/metadata/brca_annotation/metadata.txt')
metadata <- read.table(metadata_path,sep='\t',header=T)

## snpeff scoring
snpeff_scoring_path <- paste0(base_dir,'/scripts_main/hmfVariantAnnotation/data/variant_significance/snpeff_scoring.txt')
snpeff_scoring <- read.table(snpeff_scoring_path,sep='\t',header=T)

## hr genes
hr_genes <- read.table(paste0(base_dir,'/scripts_main/hmfVariantAnnotation/data/hr_genes.txt'),sep='\t',header=T)

#========= Pre-process =========#
## bind chord scores 
data <- cbind(
   sample = gene_statuses$sample,
   hrd = chord_pred[match(gene_statuses$sample,rownames(chord_pred)),'hrd'],
   response = chord_pred[match(gene_statuses$sample,rownames(chord_pred)),'response'],
   metadata[match(gene_statuses$sample, metadata$sampleId),c('cancerType','gender')],
   subset(gene_statuses,select=-sample)
); rownames(df) <- NULL

## Keep first snpeff eff (most impactful)
data$germ.snpeff_eff <- str_replace_all(data$germ.snpeff_eff,'&.+$','')
data$som.snpeff_eff <- str_replace_all(data$som.snpeff_eff,'&.+$','')

## remove TP53 rows; too many TP53 mutations
data <- data[data$gene != 'TP53',]

## remove unnecessary rows
data <- data[,!(colnames(data) %in% c('min.minor.allele.ploidy','min.copy.number','max.copy.number'))]

## add true/false positive/negative annotation
data$classif_type <- unlist(Map(function(hrd,response){
   if(hrd >= 0.5 & response != 'none'){ 'tp' }
   else if(hrd >= 0.5 & response == 'none'){ 'fp' }
   else if(hrd < 0.5 & response == 'none'){ 'tn' }
   else if(hrd < 0.5 & response != 'none'){ 'fn' }
}, data$hrd, data$response))

## simplify snpeff
insertColAfter <- function(df,pos,v,name=NULL){

   if(is.character(pos)){ pos <- which(colnames(df) == pos) }
   if(pos == 1){ stop('Just use cbind instead') }

   out <- do.call(cbind,list(
      df[,1:pos],
      v = v,
      df[,(pos+1):ncol(df)]
   ))
   
   if(!is.null(name)){ colnames(out)[pos+1] <- name }
   
   return(out)
}

toSimpleAnn <- function(v){ snpeff_scoring$ann_s2[match(v, snpeff_scoring$ann)] }

data <- insertColAfter(data,'germ.snpeff_eff',toSimpleAnn(data$germ.snpeff_eff),'germ.snpeff_eff_s')
data <- insertColAfter(data,'som.snpeff_eff',toSimpleAnn(data$som.snpeff_eff),'som.snpeff_eff_s')

subset(data, is_def==T & gene %in% c('BRCA1','BRCA2')) %>% nrow()

#========= Determine most likely cause of HRD =========#
detHrdCause <- function(){
   
   #--------- Get gene deficiencies most likely causing HRD ---------#
   ## Per sample, keep gene(s) with max hit score
   getLikelyCauses <- function(df, keep.is.def=T){
      #df = data[data$classif_type == 'fp',]
      #df = data[data$classif_type == 'tp',]
      
      df_samples <- sort(unique(df$sample))
      df_split <- lapply(df_samples, function(i){
         #i=df_samples[1]
         #i='CPCT02020740T'
         #i='CPCT02010708T'
         df_sample <- df[df$sample == i,]
         
         max_hit_score <- max(df_sample$hit_score)
         
         if(keep.is.def){
            out <- subset(df_sample, hit_score == max_hit_score | is_def)
         } else {
            out <- subset(df_sample, hit_score == max_hit_score)
         }
         
         rownames(out) <- NULL ## clean up row names
         
         return(out)
      })
      
      ## Find samples with single gene deficiency
      df_1gene <- do.call(rbind, lapply(df_split, function(i){ 
         if(nrow(i) == 1){ i } 
      })) #%>% .[order(.$hit_score, decreasing=T),]
      
      ## Find samples with multi gene deficiency (i.e. same max hit_score)
      df_multigene <- lapply(df_split, function(i){ 
         if(nrow(i) > 1){ i }
      })
      #df_multigene <- df_multigene[!sapply(df_multigene,is.null)]
      df_multigene <- do.call(rbind, df_multigene)
      
      ## Mark is_multigene_def
      df_1gene$is_multigene_def <- FALSE
      df_multigene$is_multigene_def <- TRUE
      
      df_filt <- rbind(
         df_1gene, 
         df_multigene 
      )
      
      return(df_filt)
   }
   
   data_filt <- list(
      p = getLikelyCauses( subset(data, classif_type == 'tp' | classif_type == 'fp') ),
      n = getLikelyCauses( subset(data, classif_type == 'tn' | classif_type == 'fn') )
   )
   
   #::::::::: Plots :::::::::#
   calcFreq <- function(v, name=NULL, sort=T, ...){
      #v=getHitCause(df$hit_score_origin, df$germ.snpeff_eff, df$som.snpeff_eff)
      v <- table(v, ...)
      if(sort){ v <- sort(v, decreasing=T, ...) }
      
      df <- data.frame(var = names(v), freq = as.numeric(v))
      rownames(df) <- NULL
      if(!is.null(name)){ colnames(df)[1] <- name }
      
      df$rel_freq <- df$freq/sum(df$freq)
      
      return(df)
   }
   
   #--------- Find and filter out genes are not inactivated specifically in HRD ---------#
   ## Filter out genes where (hits/total hits) in positive vs negative predicted groups is the same.
   genes_hit_filt <- lapply(data_filt, function(i){ calcFreq(i$gene, name = 'gene') })
   genes_hit_filt <- merge(genes_hit_filt$p, genes_hit_filt$n, by='gene', all=T, suffixes=c('.p','.n'))
   genes_hit_filt[is.na(genes_hit_filt)] <- 0
   
   genes_hit_filt$rel_freq_diff <- genes_hit_filt$rel_freq.p - genes_hit_filt$rel_freq.n
   
   genes_hit_filt$fisher <- (function(){
      #genes_hit_filt
      class_totals <- list(
         p = sum(genes_hit_filt$freq.p),
         n = sum(genes_hit_filt$freq.n)
      )
      
      out <- lapply(1:nrow(genes_hit_filt), function(i){
         #i=1
         r <- genes_hit_filt[i,]
         
         conting_table <- rbind(
            c(r$freq.p, class_totals$p),
            c(r$freq.n, class_totals$n)
         )
         
         fisher <- fisher.test(conting_table)
         
         return(fisher$p.value)
      })
      
      return(unlist(out))
   })()
   
   genes_hit_filt$keep <- genes_hit_filt$fisher < 0.05 & genes_hit_filt$rel_freq_diff > 0
   keep_genes <- subset(genes_hit_filt,keep==T, select=gene)[,1]
   
   #--------- genes hit x hit cause ---------#
   all_snpeff_eff <- snpeff_scoring$ann
   exist_snpeff_eff <- lapply(data_filt,function(i){ 
         i[,c('germ.snpeff_eff','som.snpeff_eff')] 
      }) %>% 
      unlist() %>% 
      unique() %>%
      .[!is.na(.)]
   
   exist_eff <- all_snpeff_eff[all_snpeff_eff %in% exist_snpeff_eff]
   exist_eff <- c('full_gene_loss', exist_eff)
   
   ## Set bar stack colors
   #col_pal <- colorRampPalette(brewer.pal(9,'Spectral'))(10)
   
   col_pal_pre <- c(
      full_gene_loss='#ff1493',
      frameshift_variant='#44146f',
      splice_acceptor_variant='#7647a2',
      splice_donor_variant='#bba3d0'
   )
   
   exist_eff_remainder <- exist_eff[!(exist_eff %in% names(col_pal_pre))]
   col_pal_gen <- colorRampPalette(brewer.pal(9,'Spectral'))(
      length(exist_eff_remainder)
   )
   col_pal <- c(col_pal_pre, col_pal_gen)
   names(col_pal) <- exist_eff
   col_pal <- rev(col_pal)
   
   plotGenesHitHitCause <- function(data, rel.freq=T, title=NULL){
      #data <- data_filt$p
      genes_hit <- calcFreq(data$gene, name='gene')

      getHitCause <- function(hit_score_origin, germ.snpeff_eff, som.snpeff_eff){
         # hit_score_origin <- data$hit_score_origin
         # germ.snpeff_eff <- data$germ.snpeff_eff
         # som.snpeff_eff <- data$som.snpeff_eff
         
         unlist(Map(function(hit_score_origin, germ.snpeff_eff, som.snpeff_eff){
            if(hit_score_origin == 'full_gene_loss'){ 
               'full_gene_loss' 
            } else if(hit_score_origin == 'loh+germ'){
               germ.snpeff_eff
            } else if(hit_score_origin == 'loh+som'){
               som.snpeff_eff
            } else if(hit_score_origin == 'germ+som'){
               c(germ.snpeff_eff,som.snpeff_eff)
            }
         }, hit_score_origin, germ.snpeff_eff, som.snpeff_eff, USE.NAMES=F))
      }
      
      genes_hit.hit_cause <- do.call(rbind,lapply(genes_hit$gene, function(i){
         #i='CFL1'
         #print(i)
         df <- data[data$gene == i,]
         cause_freq <- getHitCause(df$hit_score_origin, df$germ.snpeff_eff, df$som.snpeff_eff)
         cause_freq <- calcFreq(cause_freq,'eff')
         
         if(nrow(cause_freq) != 0){
            cause_freq <- cbind(gene = i, cause_freq)
            return(cause_freq)
         } else {
            return(NULL)
         }
         
      }))
      
      ## Re-scale relative frequency so that the top of the stacked bar is @ genes_hit$rel_freq
      genes_hit.hit_cause$rel_freq_scaled <- unlist(lapply(1:nrow(genes_hit.hit_cause),function(i){
         r <- genes_hit.hit_cause[i,]
         scale_factor <- genes_hit[genes_hit$gene == r$gene, 'rel_freq']
         scale_factor * r$rel_freq
      }))
      
      ## Keep order of gene      
      retainSortOrder <- function(v){ factor(v,levels=unique(v)) }
      genes_hit.hit_cause$gene <- retainSortOrder(genes_hit.hit_cause$gene)
      
      ## Order snpeff eff, most pathogenic first
      genes_hit.hit_cause$eff <- factor(
         genes_hit.hit_cause$eff, 
         levels=rev(c('full_gene_loss',snpeff_scoring$ann))
      )
      
      plot <- ggplot(genes_hit.hit_cause, aes(gene, if(rel.freq){rel_freq_scaled}else{freq}, fill=eff)) + 
         geom_bar(stat = "identity") +
         scale_fill_manual(values=col_pal) +
         ggtitle( ifelse(is.null(title),'Gene inactivation frequency',title) ) +
         theme(
            plot.title = element_text(hjust=0.5),
            axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
            axis.title.x = element_blank(),
            legend.title = element_blank()
         )
      
      if(rel.freq){ plot <- plot + ylab('Relative frequency') }
      else{ 
         plot <- plot + 
            ylab('Frequency') +
            scale_y_continuous(breaks= pretty_breaks())
      }
      
      return(plot)
   }
   
   #--------- Export plots ---------#
   plots <- list(
      all = plotGenesHitHitCause(rbind(data_filt$p,data_filt$n), title='Gene defiencies in whole HMF dataset'),
      
      ## Plots before fisher exact test filtering
      p = plotGenesHitHitCause(data_filt$p, title='All HR def. predicted'),
      n = plotGenesHitHitCause(data_filt$n, title='All HR prof. predicted'),
      fp = plotGenesHitHitCause(
         subset(data_filt$p,classif_type=='fp'), 
         title='False positive predictions'
      ),
      
      ## All positive samples
      filt_p = plotGenesHitHitCause(
         subset(data_filt$p, gene %in% keep_genes), 
         title='All HR def. predicted\n(fisher exact test filtered)',
         rel.freq=F
      ) + theme(legend.position = 'bottom', legend.direction='vertical'),
      
      ## False positive samples
      filt_fp = plotGenesHitHitCause(
         subset(data_filt$p, gene %in% keep_genes & classif_type=='fp'), 
         title='False positive samples\n(fisher exact test filtered)',
         rel.freq=F
      ) + theme(legend.position = 'bottom', legend.direction='vertical')
   )
   
   plot_dim <- list(
      all = c(10,5),
      p = c(10,5),
      n = c(10,5),
      fp = c(10,5),
      filt_p = c(3,6),
      filt_fp = c(2.5,5)
   )
   
   plot_out_dir <- paste0(in_dir,'/plots')
   for(i in names(plots)){
      pdf(paste0(plot_out_dir,'/gene_hit_cause_',i,'.pdf'), plot_dim[[i]][1], plot_dim[[i]][2])
      plot(plots[[i]])
      dev.off()
   }
   
   #--------- samples with multigene hits ---------#
   plotMultiHitPie <- function(){
      df <- data.frame(
         var=c('single','multi'),
         value=c(
            sum(data_filt$p$is_multigene_def == FALSE),
            subset(data_filt$p, is_multigene_def == T)$sample %>% unique() %>% length()
         )
      )
      
      plot <- ggplot(df, aes(x="", y=value, fill=var)) +
         geom_bar(width = 1, stat = "identity") +
         geom_text(label=df$value, y=cumsum(df$value)-df$value/2) +
         coord_polar("y", start=0) +
         theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid  = element_blank(),
            panel.background = element_blank(),
            legend.title = element_blank()
         )
      
      pdf(paste0(plot_out_dir,'/n_multigene_hits.pdf'), 3,3)
      plot(plot)
      dev.off()
      
   }
}

#========= Make biallelic state matrix =========#
getGeneDiplotype <- function(data, gene){
   #gene='BRCA1'
   df <- data[data$gene == gene,]
   if(nrow(df) == 0){
      out <- NULL
   } else {
      out <- do.call(rbind, lapply(1:nrow(df),function(i){
         #i=1
         r <- df[i,]
         r[is.na(r)] <- 0
         
         if(r$full_gene_loss){ 
            a1 <- a2 <- 'full_gene_loss' 
         } else if(r$loh){ 
            a1 <- 'loh'
            a2 <- 
               if(r$germ.max_score >= r$som.max_score){ r$germ.snpeff_eff_s } 
            else { r$som.snpeff_eff_s }
         } else {
            a1 <- r$germ.snpeff_eff_s
            a2 <- r$som.snpeff_eff_s
         }
         
         if(a1 == 0){ a1 <- 'no_variant' }
         if(a2 == 0){ a2 <- 'no_variant' }
         
         out <- data.frame(
            sample=r$sample,
            a1,
            a2
            # cn_break_in_gene=r$cn_break_in_gene,
            # germ.ref_loss=r$germ.ref_loss,
            # germ.alt_exists=r$germ.alt_exists,
            # som.alt_exists=r$som.alt_exists
         )
         
         return(out)
      }))
      
      colnames(out)[2:3] <- c(paste0(gene,'_1'),paste0(gene,'_2'))
   }
   return(out)
}

hr_gene_names <- 
   paste(hr_genes$string, collapse=',') %>% 
   str_remove_all(.,' ') %>% 
   str_split(.,',') %>% 
   unlist()

# l <- lapply(hr_gene_names, function(i){ 
#    message('Retrieving diplotype for: ',i)
#    getGeneDiplotype(data, i) 
# }); l <- l[!sapply(l, is.null)]

gene_diplotypes_path <- paste0(base_dir,'/scripts_main/hmfVariantAnnotation/analysis/hr_gene_def/gene_diplotypes.rds')
# saveRDS(l, gene_diplotypes_path)

#--------- Clustering ---------#
l <- readRDS(gene_diplotypes_path)
df <- Reduce(function(x,y){ merge(x,y,by='sample',all.x=T) }, l)
#rownames(df) <- df$sample; df$sample <- NULL
df[is.na(df)] <- 'no_variant' ## NA values created during merging

library(cluster)
#df <- df[1:50,2:11]
for(i in 1:ncol(df)){ df[,i] <- as.factor(df[,i]) }

dist_mat <- daisy(df[,-1])
h <- hclust(dist_mat)
h$order

#--------- Plot ---------#
df_melt <- as.data.frame(do.call(rbind, lapply(2:ncol(df), function(i){
   cbind(sample=df[[1]],allele=colnames(df)[i], eff=df[[i]])
})))

df_melt$sample <- 

## Order from least to most pathogenic
eff_order <- unique(c('no_variant',rev(snpeff_scoring$ann_s2),'loh','full_gene_loss'))
df_melt$eff <- factor(df_melt$eff, levels=eff_order[eff_order %in% df_melt$eff])

plotColors <- function(colors){
   #colors=col_pal_pre
   df <- as.data.frame(colors)
   df <- cbind(
      x=factor(rownames(df), levels=names(colors)), 
      y=1, 
      df
   ); rownames(df) <- NULL
   
   ggplot(df,aes(x,y)) + 
      geom_bar(aes(fill=x), stat='identity') +
      scale_fill_manual(values=colors) +
      theme(
         axis.title = element_blank(),
         axis.ticks = element_blank(),
         axis.text.x = element_text(angle=90,hjust=1),
         axis.text.y = element_blank(),
         panel.background = element_blank(),
         legend.position="none"
      )
}

## Predefine colors for highly pathogenic variants
col_pal_pre <- c(
   full_gene_loss='#ff1493', ## Pink
   loh='#44146f', ## Purple
   frameshift='#BD0026', ## Red
   splice_acceptor_variant='#006D2C', ## Green
   splice_donor_variant='#386CB0', ## Blue
   nonsense='#FF7F00', ## Orange
   missense='#FFFF33', ## Yellow
   'benign/likely_benign'='white',
   no_variant='white'
)
#plotColors(col_pal_pre)
###
# brewer.pal(9,'YlOrRd')
# brewer.pal(9,'Set1')
# brewer.pal(8,'Accent')
# brewer.pal(9,'BuGn')
# brewer.pal(9,'Set1')
# display.brewer.all()
###

## Generate remaining colors
# exist_eff <- levels(df_melt$eff)
# remain_eff <- exist_eff[!(exist_eff %in% names(col_pal_pre))]
# col_pal_gen <- colorRampPalette(brewer.pal(8, 'Pastel2'))( length(remain_eff) )
# names(col_pal_gen) <- remain_eff
#plotColors(col_pal_gen)

col_pal <- c(col_pal_pre,col_pal_gen)

plot <- ggplot(df_melt[1:3124*32,], aes(y=sample,x=allele)) + 
   geom_tile(aes(fill=eff)) +
   scale_fill_manual(values=col_pal_pre) +
   theme(
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
   )








