## variant_calling.R - functions for the variant calling based on bcbio RNA-Seq pipeline
##
## LSP RNAseq bcbio pipeline 
## https://github.com/chrischen1/bcbio

#' extract observed counts of different genotype from .vcf files
#'
#' @param vcf_path path contains .vcf files
#'        path must end with '/'
#' @return observed counts of different genotype from .vcf files
combine_vcf <- function(vcf_path){
  vcf_all <- NULL
  for(i in list.files(vcf_path)){
    vcf_all <- rbind(vcf_all,cbind(gsub('.vcf','',i,fixed = T),read.table(paste(vcf_path,i,sep = ''),as.is = T,header = F)))
  }
  count_info <- NULL
  for(i in 1:nrow(vcf_all)){
    trt_x <- unlist(strsplit(vcf_all$V10[i],':'))
    ctr_x <- unlist(strsplit(vcf_all$V11[i],':'))
    names(trt_x) <- names(ctr_x) <- unlist(strsplit(vcf_all$V9[i],':'))
    count_info <- rbind(count_info,c(paste(vcf_all[i,2],vcf_all[i,3],vcf_all[i,6],vcf_all[i,7],sep = ':'),as.character(vcf_all[i,1]),vcf_all[i,5],vcf_all[i,6],
                                     trt_x['RO'],trt_x['AO'],ctr_x['RO'],ctr_x['AO']))
  }
  colnames(count_info) <- c('location','sample','ref','alt','trt_count_ref','trt_count_alt','ctr_count_ref','ctr_count_alt')
  return(count_info)
}

#' filter 'PASS' variants from .vcf files
#'
#' @param vcf_path path contains .vcf files
#' @param out_path to save filtered .vcf files
#'        path must end with '/'
#' @return filtered .vcf files
vcf_filter <- function(vcf_path,out_path){
  for(i in list.files(vcf_path)){
    var_res <- NULL
    vcf  <- read.table(paste(vcf_path,i,sep = ''), quote="\"", stringsAsFactors=FALSE)
    vcf2 <- vcf[vcf$V7 == 'PASS',]
    var_res <- rbind(var_res,vcf2)
    write.table(var_res,paste(out_path,i,sep = ''),quote = F,col.names = F,row.names = F)
  }
}

#' batch coversion .vcf files to .vep files with variant effect predictor
#'
#' @param vcfDir path contains .vcf files
#' @param outDir to save  .vep files
#'        path must end with '/'
#' @return .vep files
vcf2vep <- function(vcfDir,outDir){
  bsub_fig  = 'bsub -q short -W 4:00 -o '
  vep       = 'perl ~/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl'
  dir.create(outDir,showWarnings = F)
  vcf_files <- grep('.vcf$',grep('vep',list.files(vcfDir),value = T,invert = T),value = T)
  for (i in vcf_files){
    sample_id = gsub('.vcf','',i,fixed = T)
    system(paste(bsub_fig,outDir,sample_id,'.out ',vep,' -i ',vcfDir,i,' --offline -symbol -o ',outDir,sample_id,'.vep',sep = ''))
  }
}

#'  coversion .vep files to mutation results table
#'
#' @param vep_folder path to .vep file
#' @param keep_types types of variant('CONQUENCE') that would be kept
#' @return mutation results table
vep_parser <- function(vep_folder,keep_types = c('frameshift_variant','inframe_deletion','inframe_insertion','missense_variant','protein_altering_variant','start_lost','stop_gained','stop_lost')){
  mut_mat <- NULL
  for(vep_file in list.files(vep_folder)){
    sample_id <- gsub('.vep','',vep_file,fixed = T)
    vep <- read.table(paste(vep_folder,vep_file,sep = ''), quote="\"",stringsAsFactors = F)
    vep_f <- vep[vep$V7 %in% keep_types,]
    genes <- gsub('.+SYMBOL=(.+);SYMBOL_SOURCE.+','\\1',vep_f$V14)
    mut_table <- cbind(sample_id,genes,vep_f$V4,vep_f$V5,vep_f$V7,vep_f$V2,vep_f$V3,vep_f$V8,vep_f$V9,vep_f$V10,vep_f$V11)
    mut_mat <- rbind(mut_mat,mut_table)
  }
  colnames(mut_mat) <- c('SAMPLE','SYMBOL','ENS_GENE','ENS_FEATURE','CONSEQUENCE','LOCATION','ALLELE','cDNA_POS','CDS_POS','PROT_POS','AMINO_ACID_CHANGE')
  return(mut_mat)
}