## run_de.R - Run edgeR differential expression analysis from bcbio count results and sample annotation
## Usage: Rscript run_de.R -c path/to/rnaseq.count -a path/to/annotation.txt -o path/to/output
## by Chris Chen

## Load the API
source( "de_rnaseq.R" )

## Parse command line arguments
argv <- get_args()

## argv - command line arguments 
## argv$count - Path to .count file from bcbio output
## argv$group - Path to group information file, which is a dataframe with 3 columns: group, condition and control
## argv$pairwise - If the P-values and FDR are given pairwise or as ANOVA-like test for any differences

main <- function(argv){
  pairwise <- argv$pairwise
  output <- argv$output
  counts <- read.delim(argv$count,as.is = T,row.names = 1)
  group <- read.delim(argv$annotation,as.is = T)
  if(length(grep('/$',output))==0){
    output <- paste(output,'/',sep = '')
  }
  dir.create(output,showWarnings = F)
  # remove low expression genes
  counts <- counts[apply(counts,1,min)>4,]
  de_results <- edgeR_wrapper(counts,group,!pairwise)
  if(pairwise){
    if(argv$symbol){
      cnt_sybl <- ens2symbol(rownames(de_results$logFC))
      cnt_sybl <- cnt_sybl[cnt_sybl$ensembl_gene_id!=''&cnt_sybl$hgnc_symbol!='',]
      genes <- cnt_sybl$hgnc_symbol
      names(genes) <- cnt_sybl$ensembl_gene_id
      de_results$logFC <- cbind(de_results$logFC,'symbol'=genes[rownames(de_results$logFC)])
      de_results$pmat <- cbind(de_results$pmat,'symbol'=genes[rownames(de_results$pmat)])
      de_results$fdr_mat <- cbind(de_results$fdr_mat,'symbol'=genes[rownames(de_results$fdr_mat)])
    }
    write.csv(de_results$logFC,paste(output,'logFC.csv',sep = ''))
    write.csv(de_results$pmat,paste(output,'pval.csv',sep = ''))
    write.csv(de_results$fdr_mat,paste(output,'fdr.csv',sep = ''))
  }else{
    if(argv$symbol){
      cnt_sybl <- ens2symbol(rownames(de_results))
      cnt_sybl <- cnt_sybl[cnt_sybl$ensembl_gene_id!=''&cnt_sybl$hgnc_symbol!='',]
      genes <- cnt_sybl$hgnc_symbol
      names(genes) <- cnt_sybl$ensembl_gene_id
      de_results <- cbind(de_results,'symbol'=genes[rownames(de_results)])
    }
    write.csv(de_results,paste(output,'logFC.csv',sep = ''))
  }
}

main( argv )
