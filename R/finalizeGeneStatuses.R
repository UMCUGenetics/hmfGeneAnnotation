library(magrittr)
library(stringr)
options(stringsAsFactors = F)

#========= Paths =========#
base_dir <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'
if(dir.exists('/Users/lnguyen/')){
   base_dir <- paste0('/Users/lnguyen/', base_dir)
}

#========= Load data =========#
gene_statuses <- read.table(paste0(base_dir,'/HMF_update/analysis/hr_gene_def/hmf_gene_statuses.txt'),sep='\t',header=T)

## Load CHORD output
chord_pred_path <- paste0(base_dir,'/scripts_main/mltoolkit/random_forest_training_v2/models/36_HMF_pkgMutSigs_germLt3-somEq0-lohEq0_greylistIncl_multiclass_snvContext_svContext_allSigsRel_scaleGridSearch_featPreFilt_boruta/whole_hmf_dataset_probs.txt')
#chord_pred_path <- paste0(base_dir,'/scripts_main/mltoolkit/random_forest_training_v2/models/37_HMF_pkgMutSigs_germLt3-somEq0-lohEq0_greylistIncl_multiclass_snvContext_svContext_allSigsRel_scaleGridSearch_featPreFilt_boruta_newAnn/whole_hmf_dataset_probs.txt')
chord_pred <- read.table(chord_pred_path,sep='\t',header=T)

# bind chord scores 
gene_statuses <- cbind(
   hrd = chord_pred[match(gene_statuses$sample,rownames(chord_pred)),'hrd'],
   gene_statuses
)

## Select only BRCA1/2 rows
gene_statuses_brca <- subset(gene_statuses, gene %in% c("BRCA1","BRCA2")) %>% 
   .[order(.$is_def, decreasing=T),]

#========= Compare variant annotation pipelines (Luan vs. Arne) =========#
# bind chord scores; faster
# gene_statuses_brca <- cbind(
#    hrd = unlist(lapply(gene_statuses_brca$sample, function(i){chord_pred[rownames(chord_pred) == i, 'hrd'] })),
#    gene_statuses_brca
# )

def_samples  <- list(
   luan = gene_statuses_brca %>% .[.$is_def == T,'sample'] %>% unique(),
   arne = chord_pred %>% .[.$response_simple == 1,] %>% rownames() %>% .[. != "NA"]
)

union_def_samples <- unique(unlist(def_samples))
def_samples_compare <- lapply(union_def_samples, function(i){
   #i='CPCT02010846T'
   
   data.frame(
      is_def_luan = i %in% def_samples$luan,
      is_def_arne = i %in% def_samples$arne,
      hrd = chord_pred %>% .[rownames(.) == i, 'hrd'],
      response = chord_pred %>% .[rownames(.) == i, 'response']
   )
}) %>% do.call(rbind,.)
rownames(def_samples_compare) <- union_def_samples
def_samples_compare <- def_samples_compare %>% .[order(.$hrd, decreasing=T),]

# write.table(def_samples_compare, paste0(base_dir,'/scripts_main/hmfVariantAnnotation/analysis/compare_pipelines_luan_arne/def_samples_compare.txt'),
#             sep='\t',quote=F)

#========= Finding false positives/negatives =========#
#--------- False negatives ---------#
fn_samples <- subset(def_samples_compare, is_def_luan == F & is_def_arne == T & hrd >= 0.5)

#gene_statuses_brca %>% .[.$sample %in% rownames(fn_samples),]

gene_statuses_brca %>% .[.$sample == 'CPCT02020572T',]

gene_statuses_brca %>% .[.$sample == 'CPCT02020719T',] ## BRCA1 germ mutation in original vcf; 17:41234451; was filtered out due to FILTER==SnpCluster; could also be due to nonsense in BRCA2
gene_statuses_brca %>% .[.$sample == 'CPCT02020673T',] ## BRCA2; loh_hits_ref == F; 21:25; could be due to cn_break_in_gene == T

#--------- False positives ---------#
fp_samples <- subset(def_samples_compare, is_def_luan == T & is_def_arne == F & hrd < 0.5)

subset(gene_statuses_brca, sample %in% rownames(fp_samples) & is_def == T) %>% .[order(.$hrd),]
subset(gene_statuses_brca, hit_type == 'loh+som') %>% .[order(.$hrd),]

gene_statuses_brca %>% .[.$sample == 'CPCT02040280T',] ## Correct based on snpeff; frameshift; ALT depth not taken into account
gene_statuses_brca %>% .[.$sample == 'CPCT02020680T',] ## Correct based on clinvar; Pathogenic
gene_statuses_brca %>% .[.$sample == 'CPCT02070023TII',] ## Correct based on 3 dbs; frameshift, pathogenic, pathogenic

#========= Make new BRCA response table for CHORD training =========#
def <- subset(gene_statuses_brca,is_def==T)
if(!all(table(def$sample) == 1)){ stop('Samples found with multiple gene hits') }

prof <- subset(gene_statuses_brca,is_def==F)
prof <- prof[!duplicated(prof$sample),]
prof$gene <- 'none'

annotation <- rbind(def, prof)
annotation <- annotation[!duplicated(annotation$sample),] ## rm 'none' duplicates when BRCA1 or BRCA2 is def

# write.table(
#    annotation, paste0(base_dir,'/scripts_main/hmfVariantAnnotation/analysis/hr_gene_def/hmf_brca_response.txt'),
#    sep='\t',quote=F,row.names=F
# )








