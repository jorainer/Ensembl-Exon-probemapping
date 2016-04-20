####
## getProbeSetOverlap
## function to determine the overlap of probe sets in terms of their probe composition.
## takes as input an annotation data frame with requrired columns probeset_id and probes
## the latter containing ; delimited list of probe ids.
## input paramters
## x: data frame with requrired columns probeset_id and probes
##    the latter containing ; delimited list of probe ids
getProbeSetOverlap <- function( x ){
    cat( "creating a huge matrix of ", nrow( x ), " rows by ", nrow( x ), "columns..." )
    BigMat <- matrix( nrow=nrow( x ), ncol=nrow( x ) )
    rownames( BigMat ) <- x$probeset_id
    colnames( BigMat ) <- x$probeset_id
    cat( "done\nnow start looping:\n" )
    ## i is rows, j columns
    for( i in 1:nrow( x ) ){
        i.values <- unlist( strsplit( annot[ i, "probes" ], split=";" ), use.names=FALSE )
        for( j in i:nrow( x ) ){
            BigMat[ i, j ] <- sum( i.values %in% unlist( strsplit( annot[ j, "probes" ], split=";" ), use.names=FALSE ) )
        }
        if( floor( i/500 )==ceiling(i/500) ){
            cat( "processed ", i, " of ", nrow( x ), " probe sets\n" )
        }
    }
    return( BigMat )
}
