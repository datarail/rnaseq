## Basic differential gene expression analysis via EdgeR
## Follows instructions from EdgeR Manual 1.4: Quick start
##   with minor generalization modifications
##
## LSP RNAseq bcbio pipeline 
## by Artem Sokolov, Chris Chen, et al.

## Parse command-line arguments
argv <- commandArgs( trailingOnly = TRUE )
if( length(argv) < 3 )
{
    cat( "Usage: Rscript edge.R <counts table> <metadata table> <formula>\n" )
    cat( 'Example: Rscript edge.R example/test.count example/meta.tsv "~Timepoint"\n' )
    stop( "Invalid number of arguments" )
}

cat( "Loading libraries...\n" )
suppressMessages( library( tidyverse ) )
library( stringr )
library( edgeR )

## Load all the data
cat( "Loading counts table from", argv[1], "\n" )
X <- read.delim( argv[1], row.names=1 )
cat( "Loading meta data from", argv[2], "\n" )
M <- read.delim( argv[2] )

## Identify the metadata column that best matches the counts table
v <- sapply( M, function(v) {length(intersect(v, colnames(X)))} )
if( max(v) == 0 )
    stop( "No column in the metadata file matches sample names in the counts table" )
jMatch <- which.max(v) %>% names
cat( "Using", jMatch, "column to determine a common set of samples\n" )

## Reduce to a common set of samples
M <- column_to_rownames( M, jMatch )
vCommon <- intersect( rownames(M), colnames(X) )
X <- X[,vCommon]
M <- M[vCommon,]
cat( "There are", length(vCommon), "samples in common between counts table and metadata\n" )

## Initialize the DGEList object
y <- DGEList( counts = X, samples = M )
y <- calcNormFactors(y)

## Compose the design matrix
cat( "Using the following grouping formula:", argv[3], "\n" )
mmx <- model.matrix( as.formula( argv[3] ), data = y$samples )
y <- estimateDisp( y, mmx )

## Perform quasi-likelihood F-test
qlf <- glmQLFit( y, mmx ) %>% glmQLFTest( coef=2:ncol(mmx) )

## Retrieve the results matrix
RR <- topTags( qlf, nrow(X) ) %>% as.data.frame %>% rownames_to_column( "Gene" )

## Create the differential expression directory, if needed
if( !dir.exists( "diff-expr" ) ) dir.create( "diff-expr" )

## Compose the output filename and write out the results
fn <- tools::file_path_sans_ext( argv[2] ) %>% basename %>% str_c( argv[3], ".tsv" )
fnOut <- file.path( "diff-expr", fn )
cat( "Writing results to", fnOut, "\n" )
write_tsv( RR, fnOut )
