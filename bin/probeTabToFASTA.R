probeTabToFASTA <- function( file, cols=c( "Probe.ID", "probe.x", "probe.y", "probe.sequence" ), nr=NULL ){
  cat( "reading file", file, "..." )
  chipname <- unlist( strsplit( file, split=".", fixed=TRUE ) )[ 1 ]
  Tab <- read.table( file, sep="\t", header=TRUE, as.is=TRUE )
  cat( "finished\nfile has", nrow( Tab ), "lines, ", length( unique( Tab[ , "Probe.ID" ] ) ), " unique probe ids\n" )
  cat( "generating a table with unique entries using the columns:", paste( cols, collapse=", " ), "..." )
  Tab <- unique( Tab[ , cols ] )
  ## the tab file can contain two times the header, want to exclude this
  idx <- grep( Tab[ , "probe.x" ], pattern="[a-z]", perl=TRUE )
  if( length( idx ) > 0 ){
    Tab <- Tab[ -idx,  ]
  }
  ## check length of probes...
  Lengths <-nchar( Tab[ , "probe.sequence" ] )
  idx <- which( Lengths < 25 )
  if( length( idx ) > 0 ){
    cat( "Warning! got",length(idx),"probe sequences shorter 25 nt!\n" )
    write.table( file="shorter-sequences.txt", Tab[ idx, ], sep="\t" )
    Tab <- Tab[ -idx, ]
  }
  ## check if we have probes without x and y coordinates... bad habit of
  ## Affymetrix to not have any real standardized file formats!
  noxy <- Tab[ ,"probe.x" ]==0 & Tab[ , "probe.y" ]==0
  if( any( noxy ) ){
    cat( "Warning: Got ", sum(noxy) ," probes without x and y coordinates! Trying to estimate x and y coordinates from the probe id\n" )
    if( is.null( nr ) ){
      stop( "In order to estimate x and y coordinates from the probe id I need the number of rows of the microarray to be submitted with parameter nr\n" )
    }
    require( affy )
    XY <- indices2xy( as.numeric( Tab[ noxy, "Probe.ID" ] ), nr=nr )
    Tab[ noxy, "probe.x" ] <- XY[ , "x" ]
    Tab[ noxy, "probe.y" ] <- XY[ , "y" ]
  }
  fastafile <- paste( chipname, ".fasta", sep="" )
  cat( "finished\nwriting fasta file:", fastafile, "..." )
#  con <- file( fastafile, open="w" )
  ## faster...
  Lines <- apply( Tab, MARGIN=1, function( x ){
    paste( ">probe:", chipname, ":", as.numeric( x[ "Probe.ID" ] ), ";", as.numeric( x[ "probe.x" ] ), ":", as.numeric( x[ "probe.y" ] ), "; original file:", file, "\n", x[ "probe.sequence" ], "\n", sep="" )
  })
  writeLines( Lines, con=fastafile, sep="" )
  #con <- open( File, "w" )
#  for( i in 1:nrow( Tab ) ){
#    writeLines( paste( ">probe:", chipname, ":", Tab[ i, "Probe.ID" ], ";", Tab[ i, "probe.x" ], ":", Tab[ i, "probe.y" ], "; original file:", file, sep="" ), con=con )
#    writeLines( paste( Tab[ i, "probe.sequence" ], sep="" ), con=con )
#  }
#  close( con )
  cat( "finished\nhave fun!\n" )
}
