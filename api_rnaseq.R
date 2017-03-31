library(RUVSeq)
library(edgeR)
library(biomaRt)
library(RColorBrewer)

tpm2rpkm <- function(combined,tx2gene,ercc_mapping = NULL){
  gene_mapping <- cbind('transcript'= c(tx2gene$V1,ercc_mapping$GenBank),'gene' = c(tx2gene$V2,ercc_mapping$ERCC_ID))
  genes <- gene_mapping[,2]
  names(genes) <- gene_mapping[,1]
  lib_size <- data.frame('numreads'=combined$numreads,'sample'=combined$sample)
  x <- lib_size %>% group_by(sample) %>% summarise_each(funs(sum))
  scale_factor <- x$numreads/1000000
  names(scale_factor) <- x$sample
  
  combined$RPM <- combined$numreads/scale_factor[combined$sample]
  combined$RPKM <- combined$RPM/(combined$effectiveLength/1000)
  combined$gene <- genes[combined$id]
  
  rpkm_combined <- data.frame('sample'=combined$sample,'gene'=combined$gene,'RPKM'=combined$RPKM)
  rpkm_combined_gene <- rpkm_combined %>% group_by(sample,gene)%>% summarise_each(funs(sum))
  
  rpkm_raw <- acast(rpkm_combined_gene,gene~sample)
}

# target gene processing
target_gene_processing <- function(){
  target_gene_raw <- rbind(data.frame(read.table("~/Dropbox/_ChrisProject/workspace/bcbio/annotations/control_genes.tsv", quote="\"", comment.char="", stringsAsFactors=FALSE),'source'='control',stringsAsFactors = F),
                           data.frame(read.table("~/Dropbox/_ChrisProject/workspace/bcbio/annotations/Selected_genes_201-450.tsv", quote="\"", comment.char="", stringsAsFactors=FALSE),'source'=2,stringsAsFactors = F),
                           data.frame(read.table("~/Dropbox/_ChrisProject/workspace/bcbio/annotations/Selected_genes_1-200.tsv", quote="\"", comment.char="", stringsAsFactors=FALSE),'source'=1,stringsAsFactors = F))
  x <- target_gene_raw$source
  names(x) <- target_gene_raw$V1
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  target_gene <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'hgnc_symbol', values = names(x), mart = ensembl)
  y=x[target_gene$hgnc_symbol]
  return(cbind('gene' = target_gene$ensembl_gene_id,'source'=y))
}

# get fold change, pvalue and FDR, by per cell line per time point
# order of well in samples annotation must be the same as the columns in count table
edgeR_wrapper2 <- function(cnt_time,grp_table_time,combine_fdr = F){
  design <- model.matrix(~condition,data = grp_table_time)
  y_time <- DGEList(counts=cnt_time, group=grp_table_time$condition)
  y_time <- estimateGLMCommonDisp(y_time, design)
  y_time <- estimateGLMTagwiseDisp(y_time, design)
  if(combine_fdr){
    fit <- glmFit(y_time, design)
    lrt <- glmLRT(fit, coef=2:(ncol(design)))
    lrt_tab <- topTags(lrt,n = Inf)$table[rownames(cnt_time),]
    colnames(lrt_tab) <- gsub('logFC.condition','',colnames(lrt_tab))
    return(lrt_tab)
  }
  p_mat <- fdr_mat <- logFC <- NULL
  col_names <- c()
  for(i in unique(grp_table_time$group)){
    grp_table_i <- grp_table_time[grp_table_time$group==i,]
    ctr_row <- rownames(grp_table_i)[grp_table_i$control==T]
    for (j in unique(grp_table_i$condition[grp_table_i$control!=T])){
      j_row <- rownames(grp_table_i)[grp_table_i$condition==j]
      grp_new <- rbind(grp_table_time[c(ctr_row,j_row),])
      cnt_new <- cnt_time[,rownames(grp_new)]
      result_new <- edgeR_wrapper(cnt_new,grp_new,com_disp = y_time$common.dispersion,tag_disp = y_time$tagwise.dispersion)
      if(is.null(p_mat)){
        p_mat <- result_new$pmat
        fdr_mat <- result_new$fdr_mat
        logFC <- result_new$logFC
      }else{
        p_mat <- cbind(p_mat,result_new$pmat)
        fdr_mat <- cbind(fdr_mat,result_new$fdr_mat)
        logFC <- cbind(logFC,result_new$logFC)
      }
      col_names <- c(col_names,j)
    }
  }
  colnames(p_mat) <- colnames(fdr_mat) <- colnames(logFC) <- col_names
  return(list('pmat'=p_mat,'fdr_mat'=fdr_mat,'logFC'=logFC))
}

edgeR_wrapper <- function(cnt,grp_table,w=NULL,com_disp = NULL,tag_disp = NULL){
  logFC <- NULL
  p_mat <- NULL
  fdr_mat <- NULL
  mat_label <- c()
  grp_table <- grp_table[colnames(cnt),]
  if((!is.null(w)) & (!is.null(com_disp)) & (!is.null(tag_disp))){
    design <- cbind(model.matrix(~condition,data = grp_table),w)
    y <- DGEList(counts=cnt, group=grp_table$condition)
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    com_disp <- y$common.dispersion
    tag_disp <- y$tagwise.dispersion
  }
  for (i in unique(grp_table$group)){
    selected_grp_table <- grp_table[grp_table$group==i,]
    selected_grp_table <- rbind(selected_grp_table[as.logical(selected_grp_table$control),],selected_grp_table[!as.logical(selected_grp_table$control),])
    selected_grp <- selected_grp_table$condition
    selected_cnt <- cnt[,rownames(selected_grp_table)]
    selected_w <- w[rownames(selected_grp_table),]
    y_selected <- DGEList(counts=selected_cnt, group=selected_grp)
    if(is.null(w)){
      design <- model.matrix(~condition,data = selected_grp_table)
      if((!is.null(com_disp)) & (!is.null(tag_disp))){
        y_selected$common.dispersion <- com_disp
        y_selected$tagwise.dispersion <- tag_disp
      }else{
        y_selected <- estimateGLMCommonDisp(y_selected, design)
        y_selected <- estimateGLMTagwiseDisp(y_selected, design)
      }
      y_selected <- calcNormFactors(y_selected)
      fit <- glmFit(y_selected, design)
      lrt <- glmLRT(fit, coef=2:(ncol(design)))
    }else {
      design <- cbind(model.matrix(~condition,data = selected_grp_table),selected_w)
      y_selected$common.dispersion <- com_disp
      y_selected$tagwise.dispersion <- tag_disp
      y_selected <- calcNormFactors(y_selected)
      fit <- glmFit(y_selected, design)
      lrt <- glmLRT(fit, coef=2:(ncol(design)-ncol(w)))
    }
    lrt_tab <- topTags(lrt,n = Inf)$table[rownames(cnt),]
    if(is.null(logFC)){
      logFC <- lrt_tab[,grep('logFC',colnames(lrt_tab))]
      p_mat <- lrt_tab$PValue
      fdr_mat <- lrt_tab$FDR
    }else{
      logFC <- cbind(logFC,lrt_tab[,grep('logFC',colnames(lrt_tab))])
      p_mat <- cbind(p_mat,lrt_tab$PValue)
      fdr_mat <- cbind(fdr_mat,lrt_tab$FDR)
    }
    mat_label <- c(mat_label,i)
  }
  if(!is.null(dim(p_mat))){
    colnames(p_mat) <- colnames(fdr_mat) <- mat_label
    rownames(p_mat) <- rownames(fdr_mat) <- rownames(cnt)
    colnames(logFC) <- gsub('logFC.condition','',colnames(logFC))
  }else{
    names(logFC) <- names(p_mat) <- names(fdr_mat) <- rownames(cnt)
  }
  return(list('pmat'=p_mat,'fdr_mat'=fdr_mat,'logFC'=logFC))
}

edgeR_wrapper_obs2 <- function(cnt,grp_table,w=NULL){
  logFC <- NULL
  p_mat <- NULL
  fdr_mat <- NULL
  mat_label <- c()
  grp_table <- grp_table[colnames(cnt),]
  if(is.null(w)){
    w <- matrix(0,ncol=0,nrow=ncol(cnt))
  }
  design <- cbind(model.matrix(~condition,data = grp_table),w)
  y <- DGEList(counts=cnt, group=grp_table$condition)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  com_disp <- y$common.dispersion
  tag_disp <- y$tagwise.dispersion
  y <- calcNormFactors(y)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2:(ncol(design)))
  lrt_tab_raw <- lrt$coefficients[,2:(ncol(design)-ncol(w))]
  colnames(lrt_tab_raw) <- gsub('condition','',colnames(lrt_tab_raw))
  lrt_tab <- cbind(0,lrt_tab_raw)
  colnames(lrt_tab)[1] <- grp_table$condition[1]
    
  grp_table2 <- unique(grp_table)
  cmp_table <- NULL
  for(i in unique(grp_table2$group)){
    grp_table_sub <- grp_table2[grp_table2$group==i,]
    cmp_table <- rbind(cmp_table,cbind(grp_table_sub$condition[!as.logical(grp_table_sub$control)],grp_table_sub$condition[as.logical(grp_table_sub$control)]))
  }
  for(i in 1:nrow(cmp_table)){
    if(is.null(logFC)){
      logFC <- lrt_tab[,cmp_table[i,1]]-lrt_tab[,cmp_table[i,2]]
    }else{
      logFC <- cbind(logFC,lrt_tab[,cmp_table[i,1]]-lrt_tab[,cmp_table[i,2]])
    }
  }
  colnames(logFC) <- cmp_table[,1]
  return(logFC)
}

edgeR_wrapper_obs <- function(cnt,samples,time =c(6,24),w=NULL,com_disp = NULL,tag_disp=NULL){
  grp <- paste(samples$CellLine,samples$DrugName,samples$Conc,samples$Time,sep = '_')
  logFC <- NULL
  p_mat <- NULL
  fdr_mat <- NULL
  mat_label <- c()
  if(is.null(w)){
    design <- model.matrix(~grp)
  }else {
    design <- cbind(model.matrix(~grp),w)
  }
  # get overall disp and use them for per cell line
  y <- DGEList(counts=cnt, group=grp)
  if(is.null(com_disp)){
    y <- estimateGLMCommonDisp(y, design)
    com_disp <- y$common.dispersion
  }
  if(is.null(tag_disp)){
    y <- estimateGLMTagwiseDisp(y, design)
    tag_disp <- y$tagwise.dispersion
  }
  for (i in unique(samples$CellLine)){
    for (j in time){
      selected_col <- samples$CellLine==i & (samples$Time == 0 | samples$Time == j)
      selected_cnt <- cnt[,selected_col]
      selected_grp <- grp[selected_col]
      y_selected <- DGEList(counts=selected_cnt, group=selected_grp)
      y_selected$common.dispersion <- com_disp
      y_selected$tagwise.dispersion <- tag_disp
      if(is.null(w)){
        design_selected <- model.matrix(~selected_grp)
        fit <- glmFit(y_selected, design_selected)
        lrt <- glmLRT(fit, coef=2:ncol(design_selected))
      }else{
        selected_w <- data.frame(w)[selected_col,]
        design_selected <- cbind(model.matrix(~selected_grp),selected_w)
        
        fit <- glmFit(y_selected, design_selected)
        lrt <- glmLRT(fit, coef=2:(ncol(design_selected)-1))
      }
      lrt_tab <- topTags(lrt,n = Inf)$table[rownames(cnt),]
      if(is.null(logFC)){
        logFC <- lrt_tab[,grep('logFC',colnames(lrt_tab))]
        p_mat <- lrt_tab$PValue
        fdr_mat <- lrt_tab$FDR
      }else{
        logFC <- cbind(logFC,lrt_tab[,grep('logFC',colnames(lrt_tab))])
        p_mat <- cbind(p_mat,lrt_tab$PValue)
        fdr_mat <- cbind(fdr_mat,lrt_tab$FDR)
      }
      mat_label <- c(mat_label,paste(i,j,sep = '_'))
    }
  }
  colnames(p_mat) <- colnames(fdr_mat) <- mat_label
  rownames(p_mat) <- rownames(fdr_mat) <- rownames(cnt)
  colnames(logFC) <- gsub('logFC.selected_grp','',colnames(logFC))
  return(list('pmat'=p_mat,'fdr_mat'=fdr_mat,'logFC'=logFC))
}

ens2symbol <- function(ens){
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  target_gene <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id', values = ens, mart = ensembl)
  return(target_gene)
}


get_fc_exact <- function(cnt,grp,x='logFC'){
  DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
  DGEx1f <- calcNormFactors(DGEx1)
  design <- model.matrix(~grp)
  DGEx1f <- estimateDisp(DGEx1f,design)
  logFC <- NULL
  col_label <- c()
  for(i in unique(samples$CellLine)){
    grp_cell <- grep(i,grp,value = T)
    grp_cell_ctr <- intersect(grep('-',grp,value = T),grp_cell)
    grp_cell_trt <- grp_cell[grp_cell!=grp_cell_ctr]
    for(j in grp_cell_trt){
      col_label <- c(col_label,j)
      fit <-exactTest(DGEx1f, pair=c(grp_cell_ctr,j)) 
      tags <- topTags(fit, n=Inf)$table
      val <- which(x==colnames(tags))
      tags2 <- tags[,val]
      names(tags2) <- tags$genes
      tags2 <- tags2[DGEx1$genes$genes]
      if(is.null(logFC)) logFC <- tags2
      else logFC <- cbind(logFC,tags2)
    }
  }
  colnames(logFC) <- col_label
  return(logFC)
}

get_fc_glm <- function(cnt,grp){
  DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
  DGEx1f <- calcNormFactors(DGEx1)
  logFC <- NULL
  for(i in unique(samples$CellLine)){
    grp_cell <- grep(i,grp)
    grp_cell_ctr <- intersect(grep('-',grp),grp_cell)
    grp_cell_trt <- grp_cell[grp_cell!=grp_cell_ctr]
    grp_cell_all <- c(grp_cell_ctr,grp_cell_trt)
    y_d = DGEx1f[,grp_cell_all]
    grp_c <- y_d$samples$group
    design <- model.matrix(~grp_c)
    y_d <- estimateDisp(y_d, design, robust=TRUE)
    fit <- glmFit(y_d, design)
    lrt <- glmLRT(fit, coef=2:ncol(fit$design))
    tags <- topTags(lrt, n=Inf)$table
    tags2 <- tags[,2:ncol(fit$design)]
    rownames(tags2) <- tags$genes
    tags2 <- tags2[DGEx1$genes$genes,]
    if(is.null(logFC)) logFC <- tags2
    else logFC <- cbind(logFC,tags2)
  }
  colnames(logFC) <- gsub('logFC.grp_c','',colnames(logFC))
  return(logFC)
}

get_fc_qglm <- function(cnt,grp){
  DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
  DGEx1f <- calcNormFactors(DGEx1)
  logFC <- NULL
  for(i in unique(samples$CellLine)){
    grp_cell <- grep(i,grp)
    grp_cell_ctr <- intersect(grep('-',grp),grp_cell)
    grp_cell_trt <- grp_cell[grp_cell!=grp_cell_ctr]
    grp_cell_all <- c(grp_cell_ctr,grp_cell_trt)
    y_d = DGEx1f[,grp_cell_all]
    grp_c <- y_d$samples$group
    design <- model.matrix(~grp_c)
    y_d <- estimateDisp(y_d, design, robust=TRUE)
    fit <- glmQLFit(y_d, design)
    lrt <- glmQLFTest(fit, coef=2:ncol(fit$design))
    tags <- topTags(lrt, n=Inf)$table
    tags2 <- tags[,2:ncol(fit$design)]
    rownames(tags2) <- tags$genes
    tags2 <- tags2[DGEx1$genes$genes,]
    if(is.null(logFC)) logFC <- tags2
    else logFC <- cbind(logFC,tags2)
  }
  colnames(logFC) <- gsub('logFC.grp_c','',colnames(logFC))
  return(logFC)
}

get_raw_fc <- function(cnt,samples){
  cnt <- cnt+0.01
  raw_fc <- NULL
  for (i in unique(samples$CellLine)){
    well_ctr <- samples$well[samples$CellLine==i & samples$ctrl]
    well_trt <- samples$well[samples$CellLine==i & (!samples$ctrl)]
    if(is.null(raw_fc)) raw_fc <- cnt[,well_trt]/cnt[,well_ctr]
    else raw_fc <- cbind(raw_fc,cnt[,well_trt]/cnt[,well_ctr])
  }
  return(log(raw_fc))
}

getFC_ruv <- function(cnt,grp,w = NULL){
  DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
  DGEx1f <- calcNormFactors(DGEx1)
  if(is.null(w)) design <- model.matrix(~grp)
  else design <- model.matrix(~grp+w)
  y_d <- estimateDisp(DGEx1f, design, robust=TRUE)
  fit <- glmFit(y_d, design)
  lrt <- glmLRT(fit, coef=2:ncol(fit$design))
  tags <- topTags(lrt, n=Inf)$table
  lfc_ind <- grep('logFC',colnames(tags))
  lfc_ind <- lfc_ind[colnames(tags)[lfc_ind] != 'logFC.w']
  tags2 <- tags[,lfc_ind]
  rownames(tags2) <- tags$genes
  tags2 <- tags2[rownames(cnt),]
  colnames(tags2) <- gsub('logFC.grp','',colnames(tags2))
  tags3 <- cbind('_._'=0,tags2)
  ctr_ind <- grep('_._',colnames(tags3),fixed = T)
  FC_table <- NULL
  for(i in ctr_ind){
    FC_table <- rbind(FC_table,cbind(rep(i,6),(i+1):(i+6)))
  }
  logFC <- NULL
  for (i in 1:nrow(FC_table)){
    fc_diff <- tags3[,FC_table[i,2]]-tags3[,FC_table[i,1]]
    if(is.null(logFC)) logFC <- fc_diff
    else logFC <- cbind(logFC,fc_diff)
  }
  colnames(logFC) <- colnames(tags3)[-ctr_ind]
  rownames(logFC) <- rownames(tags3)
  return(logFC) 
}


#p-values
get_pval_ruv <- function(cnt,grp,w = NULL){
  DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
  DGEx1f <- calcNormFactors(DGEx1)
  if(is.null(w)) {
    design <- model.matrix(~grp)
  }else {
    design <- model.matrix(~grp+w)
  }
  y_d <- estimateDisp(DGEx1f, design, robust=TRUE)
  fit <- glmFit(y_d, design)
  pvals <- NULL
  for(i in 2:ncol(fit$design)){
    lrt <- glmLRT(fit, coef=i)
    tags <- topTags(lrt, n=Inf)$table
    tags2 <- tags$PValue
    names(tags2) <- tags$genes
    if(is.null(pvals)){
      pvals <- tags2
    }else{
      pvals <- cbind(pvals,tags2)
    }
  }
  colnames(pvals) <- colnames(fit$design)[-1]
  return(pvals) 
}

getFC <- function(cnt,grp,hyper_grp=NULL,w = NULL){
  DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
  DGEx1f <- calcNormFactors(DGEx1)
  if(is.null(hyper_grp)){
    hyper_grp <- rep(1,length(grp))
  }else{
    hyper_grp <- hyper_grp
  }
  logFC <- NULL
  for(i in unique(hyper_grp)){
    DGEx1f_sub <- DGEx1f[,hyper_grp ==i]
    grp_sub <- grp[hyper_grp ==i]
    if(is.null(w)) {
      design <- model.matrix(~grp_sub)
      y_d <- estimateDisp(DGEx1f_sub, design, robust=TRUE)
    }else{
      w_sub <- w[hyper_grp ==i]
      design0 <- model.matrix(~grp_sub)
      design <- model.matrix(~grp_sub+w_sub)
      y_d <- estimateDisp(DGEx1f_sub, design0, robust=TRUE)
    }
    fit <- glmFit(y_d, design)
    lrt <- glmLRT(fit, coef=2:length(unique(grp_sub)))
    tags <- topTags(lrt, n=Inf)$table
    tags2 <- tags[,2:length(unique(grp_sub))]
    rownames(tags2) <- tags[,1]
    colnames(tags2) <- gsub('logFC.grp_sub','',colnames(tags2))
    if(is.null(logFC)){
      logFC <- tags2[rownames(cnt),]
    }else{
      logFC <- cbind(logFC,tags2[rownames(cnt),])
    }
  }
  return(logFC)
}
