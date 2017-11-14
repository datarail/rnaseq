## Merges individual Salmon output files
##
## LSP RNAseq bcbio pipeline 
## by Artem Sokolov, Chris Chen, et al.

cat( "Loading libraries...\n" )
library( stringr )

suppressMessages( library( tidyverse ) )
suppressMessages( library( magrittr ) )

## Identify all salmon files
sf.files <- list.files( pattern="quant\\.sf$", recursive=TRUE )
sf.files <- sf.files[ grep( "final", sf.files ) ]
cat( "Found the following salmon files:", sf.files, sep="\n" )
cat( "-------------------------------\n" )

## Infer sample names from the identified filenames
sspl <- str_split( sf.files, "/" )
jns <- which( sapply( sspl, nth, -2 ) != "salmon" )
if( length(jns) > 0 )
{
    cat( "Warning: the following quant.sf files are not in a salmon directory.\n" )
    cat( "This can lead to improper sample name inference\n" )
    cat( "Consider excluding them from the merger\n" )
    cat( sf.files[jns], sep="\n" )
    cat( "-------------------------------\n" )
}
sampleIDs = sapply( sspl, nth, -3 )

## Infer the longest common path
jCommon <- sapply( sspl, length) %>% min %>% seq( 1, . ) %>%
    lapply( sspl, "[", . ) %>% do.call( cbind, . ) %>%
    apply( 1, function(v) {length(unique(v))} ) %>%
    is_greater_than( 1 ) %>% which %>% min %>% subtract( 1 )
if( jCommon == 0 ) {
    dirOut <- "."
} else {
    dirOut <- sspl[[1]][1:jCommon] %>% as.list %>% do.call( file.path, . )
}
fnCombined <- file.path( dirOut, "combined.sf" )

## Load individual files and combine into a common matrix
XX <- lapply( 1:length(sf.files), function(i) {
    cat( "Loading", sf.files[i], "\n" )
    suppressMessages(read_tsv( sf.files[i] )) %>%
        mutate( Sample = sampleIDs[i] ) }) %>%
    bind_rows

## Store the combined output
cat( "Writing", fnCombined, "\n" )
write_tsv( XX, fnCombined )

## Identify and load transcript -> gene name mapping
fn.t2g <- list.files( pattern="tx2gene\\.csv$", recursive=TRUE )
stopifnot( length(fn.t2g) == 1 )
T2G <- suppressMessages(read_csv( fn.t2g, col_names=FALSE )) %>%
    rename( Name = X1, Gene = X2 )

## Map transcripts to genes and compute a summary statistic for each gene, sample pair
## The summary statistic is taken to be the sum across all transcripts belonging to the same gene
cat( "Summarizing counts and TPM values across each gene...\n" )
SS <- inner_join( T2G, XX, by="Name" ) %>% select( Gene, Sample, TPM, NumReads ) %>%
    group_by( Gene, Sample ) %>% summarize_all( funs(sum) ) %>% ungroup

## Compose the output filenames
fnCounts <- file.path( dirOut, "sf_combined_counts.tsv" )
fnTPM <- file.path( dirOut, "sf_combined_tpm.tsv" )

## Write out the counts and TPM tables
cat( "Writing", fnCounts, "\n" )
SS %>% select( Gene, Sample, NumReads ) %>% spread( Sample, NumReads ) %>% write_tsv( fnCounts )
cat( "Writing", fnTPM, "\n" )
SS %>% select( Gene, Sample, TPM ) %>% spread( Sample, TPM ) %>% write_tsv( fnTPM )

