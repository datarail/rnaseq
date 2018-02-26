## Basic differential gene expression analysis via EdgeR
## Follows instructions from EdgeR Manual 1.4: Quick start
##   with minor generalization modifications
##
## LSP RNAseq bcbio pipeline 
## by Artem Sokolov, Chris Chen, et al.

## Parse command-line arguments
argv <- commandArgs( trailingOnly = TRUE )
if( length(argv) < 4 )
{
    cat( "Usage: Rscript edge.R <counts table> <metadata table> <contrast column> <control value>\n" )
    cat( 'Example: Rscript edge.R example/test.count example/meta.tsv Timepoint 0h\n' )
    stop( "Invalid number of arguments" )
}

cat( "Loading libraries...\n" )
suppressMessages( library( tidyverse ) )
suppressMessages( library( edgeR ) )
library( stringr )

## Load all the data
cat( "Loading counts table from", argv[1], "\n" )
X <- read.delim( argv[1], row.names=1, check.names=FALSE )
cat( "Loading meta data from", argv[2], "\n" )
M <- read.delim( argv[2], as.is = TRUE )

## Argument verification
if( !(argv[3] %in% colnames(M)) )
    stop( "No such column ", argv[3], " in the meta file" )
if( !(argv[4] %in% M[[argv[3]]] ) )
    stop( "Control value ", argv[4], " never occurs in column ", argv[3] )

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

## Convert the desired grouping column to factor
## Use control as the first level in this factor
cat( "Using", argv[4], "as control value for contrasts along column", argv[3], "\n" )
vGrp <- M[[argv[3]]]
lvl <- c(argv[4], vGrp) %>% unique
M <- M %>% mutate( UQ(argv[3]) := factor( vGrp, lvl ) )

## Initialize the DGEList object
y <- DGEList( counts = X, samples = M )
y <- calcNormFactors(y)

## Compose the design matrix
f <- str_c( "~", argv[3] )
cat( "Using the following formula for setting up contrasts:", f, "\n" )
mmx <- model.matrix( as.formula( f ), data = y$samples )
y <- estimateDisp( y, mmx )

## Perform quasi-likelihood F-test
qlf <- glmQLFit( y, mmx ) %>% glmQLFTest( coef=2:ncol(mmx) )

## Retrieve the results matrix
RR <- topTags( qlf, nrow(X) ) %>% as.data.frame %>% rownames_to_column( "Gene" )

## Create the differential expression directory, if needed
if( !dir.exists( "diff-expr" ) ) dir.create( "diff-expr" )

## Compose the output filename and write out the results
fn <- tools::file_path_sans_ext( argv[2] ) %>% basename %>% str_c( "-", argv[3], ".tsv" )
fnOut <- file.path( "diff-expr", fn )
cat( "Writing results to", fnOut, "\n" )
write_tsv( RR, fnOut )
