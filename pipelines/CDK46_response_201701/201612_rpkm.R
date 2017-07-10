source('../common/de_rnaseq.R')
# change working directory to  201612_de/final/runDate_201612_de/
sf_file = 'combined.sf'
tx_file = 'tx2gene.csv'

sf <- read.delim(sf_file,as.is = T)
tx <- read.csv(tx_file,as.is = T,header=FALSE)
sample_info <- read.delim('../../DATA/RNA_CDK46_perturbation/201612_annotation.tsv',as.is = T)

samples <- paste(sample_info$Cell.line,sample_info$Treatment,sample_info$Time.point)
names(samples) <- sample_info$Sample
samples <- gsub(' N/A N/A','',samples)
rpkm <- tpm2rpkm(sf,tx)
colnames(rpkm) <- samples[gsub('.+S(\\d+)_.+','\\1',colnames(rpkm))]

genes <- ens2symbol(rownames(rpkm))
genes <- genes[genes$hgnc_symbol!='' & !(genes$ensembl_gene_id %in% genes$ensembl_gene_id[duplicated(genes$ensembl_gene_id)]) & 
                 !(genes$hgnc_symbol %in% genes$hgnc_symbol[duplicated(genes$hgnc_symbol)]),]
gene_name <- genes$hgnc_symbol
names(gene_name) <- genes$ensembl_gene_id
rpkm_symbol <- rpkm[rownames(rpkm)%in%names(gene_name),]
rownames(rpkm_symbol) <- gene_name[rownames(rpkm_symbol)]

write.table(rpkm,'../../DATA/RNA_CDK46_perturbation/201612_rpkm_ens.tsv',sep = '\t')
write.table(rpkm_symbol,'../../DATA/RNA_CDK46_perturbation/201612_rpkm_symbol.tsv',sep = '\t')
