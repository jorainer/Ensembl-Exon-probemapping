## Compare CDFs a and b with each other.
compareCDF <- function( a, b ){
    ## what do we want to compare:
    ## No of probe sets.
    ## No of unique probes.
    ## Average number of probes per probe set.
    ## % of probe sets with identical probe content.
    require( affxparser )
    A <- readCdfUnits( a )
    B <- readCdfUnits( b )
    cat( "-----------------------------------------\n> Number of defined probe sets:\n" )
    cat( " ", a, ": ", length( A ),"\n" )
    cat( " ", b, ": ", length( B ),"\n" )
    ##
    cat( "\n> Number of unique probes:\n" )
    Aprobes <- lapply( A, function( x ){
        return( paste( x$groups[[1]]$x, x$groups[[1]]$y, sep=":" ) )
    } )
    cat( " ", a, ": " )
    cat( length( unique( unlist( Aprobes ) ) ), "\n" )
    Bprobes <- lapply( B, function( x ){
        return( paste( x$groups[[1]]$x, x$groups[[1]]$y, sep=":" ) )
    } )
    cat( " ", b, ": " )
    cat( length( unique( unlist( Bprobes ) ) ), "\n" )
    ##
    cat( "\n> Average number of probes per probe set:\n" )
    cat( " ", a, ": " )
    cat( mean( unlist( lapply( Aprobes, length ) ) ), "\n" )
    cat( " ", b, ": " )
    cat( mean( unlist( lapply( Bprobes, length ) ) ), "\n" )
    ##
    cat( "\n> Proportion of probe sets with identical probe composition:\n" )
    Aprobes <- unlist( lapply( Aprobes, function( x ){
        paste( sort(x), collapse=";" )
    } ) )
    Bprobes <- unlist( lapply( Bprobes, function( x ){
        paste( sort(x), collapse=";" )
    } ) )
    cat( length( intersect( Aprobes, Bprobes ) ), "\n" )
    ## get the list of probe sets with different probe composition.
    return( Aprobes[ !( Aprobes %in% Bprobes ) ] )
}
