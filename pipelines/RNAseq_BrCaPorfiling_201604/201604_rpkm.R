source('../common/de_rnaseq.R')
#get combined.sf if not exist
# change working directory to  201604_de/work/
# sf_files <- list.files(path="./",pattern='*\\.sf',recursive=T)
# sf_info <- NULL
# for(i in sf_files){
#   si <- read.delim(i,as.is = T)
#   si <- cbind(si,'sample'=gsub('salmon/(.+)/quant/.+','\\1',i),'id'=si$Name)
#   sf_info <- rbind(sf_info,si)
# }
# sf <- sf_info
# colnames(sf) <- c('name','length','effectiveLength','tpm','numreads','sample','id')
# change this according to run date
# write.table(sf,'../final/runDate_201604_de/combined.sf',sep='\t')

# change working directory to  201604_de/final/runDate_201604_de/
sf_file = 'combined.sf'
tx_file = 'tx2gene.csv'

sf <- read.delim(sf_file,as.is = T)
tx <- read.csv(tx_file,as.is = T,header=FALSE)
rpkm <- tpm2rpkm(sf,tx)

sample_info <- read.delim('../../DATA/RNA_baseline/201604_annotation.tsv',as.is = T)
samples <- sample_info$cellines
rep_tag <- sample_info$rep_tag
names(samples) <- names(rep_tag) <- sample_info$sample_id
colnames(rpkm) <- gsub('.bio.replic','',colnames(rpkm),fixed = T)
colnames(rpkm) <- gsub('.tech.replic','',colnames(rpkm),fixed = T)
rep_tag_insert <- rep_tag[gsub('X','',colnames(rpkm))]
colnames(rpkm) <- samples[gsub('X','',colnames(rpkm))]

genes <- ens2symbol(rownames(rpkm))
genes <- genes[genes$hgnc_symbol!='' & !(genes$ensembl_gene_id %in% genes$ensembl_gene_id[duplicated(genes$ensembl_gene_id)]) & 
                 !(genes$hgnc_symbol %in% genes$hgnc_symbol[duplicated(genes$hgnc_symbol)]),]
gene_name <- genes$hgnc_symbol
names(gene_name) <- genes$ensembl_gene_id
rpkm_symbol <- rpkm[rownames(rpkm)%in%names(gene_name),]
rownames(rpkm_symbol) <- gene_name[rownames(rpkm_symbol)]

rpkm_ens <- rbind(rpkm,rep_tag_insert)
rpkm_ens <- rbind(rpkm_ens[nrow(rpkm_ens),],rpkm_ens[-nrow(rpkm_ens),])
rownames(rpkm_ens)[1] <- 'rep_tag'

rpkm_symbol <- rbind(rpkm_symbol,rep_tag_insert)
rpkm_symbol <- rbind(rpkm_symbol[nrow(rpkm_symbol),],rpkm_symbol[-nrow(rpkm_symbol),])
rownames(rpkm_symbol)[1] <- 'rep_tag'

write.table(rpkm_symbol,'../../DATA/RNA_baseline/201604_rpkm_symbol.tsv',sep = '\t')

#merge replicates
rpkm_merge_rep <- NULL
sample_names <- c()
for(i in unique(colnames(rpkm))){
  sample_names <- c(sample_names,i)
  if(sum(colnames(rpkm)==i)==1){
    rpkm_merge_rep <- rbind(rpkm_merge_rep,rpkm[,colnames(rpkm)==i])
  }else{
    rpkm_merge_rep <- rbind(rpkm_merge_rep,apply(rpkm[,colnames(rpkm)==i],1,mean))
  }
}
rownames(rpkm_merge_rep) <- sample_names
rpkm_merge_rep <- t(rpkm_merge_rep)
rpkm_merge_rep_symbol <- rpkm_merge_rep[rownames(rpkm_merge_rep)%in%names(gene_name),]
rownames(rpkm_merge_rep_symbol) <- gene_name[rownames(rpkm_merge_rep_symbol)]

write.table(rpkm_merge_rep,'../../DATA/RNA_baseline/201604_rpkm_merge_rep_ens.tsv',sep = '\t')
write.table(rpkm_merge_rep_symbol,'../../DATA/RNA_baseline/201604_rpkm_merge_rep_symbol.tsv',sep = '\t')
