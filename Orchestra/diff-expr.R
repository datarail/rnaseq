## diff-expr.R - main wrapper for the differential expression analysis
##
## LSP RNAseq pipeline
## https://github.com/sorgerlab/rnaseq

#' Main wrapper for the differential expression analysis
#'
#' Applies DESeq2 or edgeR to the provided count table and metadata along the attribute of interest.
#'
#' @param X p-by-n count matrix for p genes across n samples.
#' @param Y n-by-k matrix of k attributes for the n samples. rownames(Y) must match colnames(X)
#' @param dattr string defining the differential attribute of interest
#' @param tattr (option) string defining the temporal attribute
#' @param method one of {"DESeq2", "edgeR"} defining which differential expression method to apply
#' @return Matrix of dimension p-by-... that contains differential expression for each gene
diff.expr <- function( X, Y, dattr, tattr = NULL, method = "DESeq2" )
{
    ## ######################
    ## Parameter verification
    ## ######################
    if( all( rownames(Y) == colnames(X) ) == FALSE )
        stop( "Row names of Y must match column names of X." )

    if( dattr %in% colnames(Y) == FALSE )
        stop( "dattr must be an attribute of Y." )

    if( is.null( tattr ) == FALSE && tattr %in% colnames(Y) == FALSE )
        stop( "tattr must be an attribute of Y or NULL." )

    ## #############################################################
    ## Passing the call to the corresponding DESeq2 / edgeR function
    ## #############################################################

    ## Additional notes:
    ## - If tattr is NULL, DESeq2 / edgeR should run standard differential pipeline.
    ##   Otherwise, use the following guides to construct pairwise comparison across all timepoints:
    ##   DESeq2: http://www.bioconductor.org/help/workflows/rnaseqGene/#time-course-experiments
    ##   edgeR: Section 3.3.4 in https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    ## - If dattr is multinomial, set up a one-vs-one or one-vs-all scheme
    
    if( method == "DESeq2" )
    {
        ## Call DESe2 wrapper here
    } else if ( method == "edgeR" ) {
        ## Call edgeR wrapper here
    } else
        stop( "Unrecognized method. Must be \"DESeq2\" or \"edgeR\"." )
}
