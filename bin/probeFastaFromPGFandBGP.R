## This function extracts probe sequence, id and x and y positions on the array from .pgf and .bgp files from
## the Library zip files from Affymetrix.
## Ideally, we would wish to get this information directly from Affymetrix, but their probe sequence files contain much less
## probes...
## Note that we are "estimating" the x and y position from the probe id and the number of rows of the array.
## nr shold be the number of rows of the array; has only to be specified if it can not be found within the file (as usual Affy
## does not have standardized file formats...).
## nr for HuGene, MoGene 1.0: 1050
## nr for HuGene 2.0: 1612
## nr for HuEx, MoEx 1.0: 2560
## CAUTION: probe sequences in the pgf and bgp are 3' to 5' instead of 5' - 3'! see http://media.affymetrix.com/support/developer/powertools/changelog/file-format-pgf.html
probeFastaFromPGFandBGP <- function( pgf, nr, chiptype="st" ){
    require( affy )
    require( Biostrings )
    if( missing( pgf ) ){
        stop( "The pgf file has to be defined!" )
    }
    if( chiptype=="st" ){
        seqtransform <- function( x ){
            return( as.character( reverse( DNAString( x ) ) ) )
        }
    }else{
        seqtransform <- function( x ){
            return( as.character( reverse( DNAString( x ) ) ) )
        }
    }
    con <- file( pgf, open="r" )
    SequenceTable <- matrix( ncol=4, nrow=1 )
    colnames( SequenceTable ) <- c( "probe_id", "x", "y", "sequence" )
    ## have to define these
    outfile <- ""
    outcon <- ""
    chip_type <- ""
    nr_rows <- 0
    nr_cols <- 0
    create_date <- ""
    ## get the file name in case a path was submitted...
    pgf_file <- unlist( strsplit( pgf, "/" ) )
    pgf_file <- pgf_file[ length( pgf_file) ]
    counter <- 0
    ##
    tempfile <- "temp.out"
    unlink( tempfile )
    outcon <- file( tempfile, open="a" )
    ## read the file line by line and add to a data.frame. that way we can then check for probe sequence duplications..
    cat( "processing file: ", pgf_file, "..." )
    while( length( aLine <- readLines( con, n=1, warn=FALSE ) ) > 0 ){
        if( length( grep( aLine, pattern="^#" ) ) > 0 ){
            ## got the header.
            Vals <- unlist( strsplit( aLine, split="=" ) )
            if( Vals[ 1 ]=="#%chip_type" ){
                chip_type <- Vals[ 2 ]
                outfile <- paste0( chip_type, ".fasta" )
            }
            if( Vals[ 1 ]=="#%create_date" ){
                create_date <- Vals[ 2 ]
            }
            if( Vals[ 1 ]=="#%num-cols" ){
                nr_cols<- as.numeric( Vals[ 2 ] )
            }
            if( Vals[ 1 ]=="#%num-rows" ){
                nr_rows<- as.numeric( Vals[ 2 ] )
                if( nr_rows!=nr_cols) stop( paste0( "Argh, got ", nr_rows, " rows but ", nr_cols, " cols, don't know what to do!" ))
            }
            if( Vals[ 1 ]=="#%header0" ){
                ## if we don't have nr rows/nr cols here, stop!
                if( nr_rows==0 & missing( nr ) ){
                    close( con )
                    close( outcon )
                    stop( "No information about number of probes found in the file! Number of rows has to be manually specified with the nr parameter!" )
                }
                if( nr_rows==0 ){
                    nr_rows <- nr
                }
            }
        }else{
            ## read the line and append probe sequence, x and y to outfile
            ## x and y will be estimated on the probe id.
            ## we are looking for lines that start with two empty entries or have 8 entries.
            ## these are then: "", "", probe_id, type, gc_count, probe_length, interrogation_position, probe_sequence
            ## we will only extract probe_id and sequence from that...
            Vals <- unlist( strsplit( aLine, split="\t" ) )
            if( length( Vals )==8 ){
                ## if probe length!=25, skip.
                if( as.numeric( Vals[ 6 ] )==25 ){
                    ## probe id: 3
                    ## sequence: 8
                    xy <- indices2xy( as.numeric( Vals[ 3 ] ), nr=nr_rows)
                    #SequenceTable <- rbind( SequenceTable, c( Vals[ 3 ], xy, Vals[ 8 ] ) )
                    writeLines( paste0( Vals[ 3 ], "\t", xy[ 1 ], "\t", xy[ 2 ], "\t", Vals[ 8 ] ), con=outcon )
                    #writeLines( c( paste0( ">probe:", chip_type, ":", Vals[ 3 ], ";", xy[ 1 ], ":", xy[ 2 ], "; original file:", pgf_file, ";create date ", create_date ), Vals[ 8 ] ), con=outcon )
                }else{
                    warning( "Got a probe with length ", Vals[ 6 ], ": ", Vals[ 8 ] )
                }
            }
        }
    }
    close( con )
    ## bg file(s)
    Path <- unlist( strsplit( pgf, split="/" ) )
    if( length( Path ) > 1 ){
        Path <- paste( Path[ 1:(length( Path )-1) ], collapse="/" )
    }else{
        Path <- "."
    }
    Files <- dir( Path, pattern=sub( pgf_file, pattern=".pgf", replacement="", fixed=TRUE ) )
    bg_files <- Files[ grep( Files, pattern=".bgp$" ) ]
    for( bg_file in bg_files ){
        cat("done\nprocessing background probe file: ", bg_file, "...")
        ## now getting probes from the background probe file.
        ## read the background probes and append to file.
        BGSeq <- read.table( paste( c( Path, bg_file ), collapse="/" ), header=TRUE, as.is=TRUE, sep="\t" )
        for( i in 1:nrow( BGSeq ) ){
            if( as.numeric( BGSeq[ i, "probe_length" ] ) == 25 ){
                writeLines( paste( BGSeq[ i, c( "probe_id", "x", "y", "probe_sequence" ) ], collapse="\t" ), con=outcon )
            }
        }
    }
    close( outcon )  ## that's the temp file here.
    SequenceTable <- read.table( tempfile, sep="\t", as.is=TRUE )
    cat( "done\nchecking probe sequences for duplicates etc:" )
    nr_orig <- nrow( SequenceTable )
    SequenceTable <- unique( SequenceTable )
    cat( "Had originally ", nr_orig, " sequences, after removing duplicates (based on id, x, y, sequence): ", nrow( SequenceTable ),"\n")
    cat( "writing output: ", outfile, "..." )
    ## writing fasta
    unlink( outfile )
    outcon <- file( outfile, open="a" )
    for( i in 1:nrow( SequenceTable ) ){
        writeLines( paste0( ">probe:", chip_type, ":", SequenceTable[ i, 1 ], ";", SequenceTable[ i, 2 ], ":", SequenceTable[ i, 3 ], "; original file:", pgf_file, ";create date ", create_date, "\n", seqtransform( SequenceTable[ i, 4 ] ), "\n" ), sep="", con=outcon )
        if( floor(i/500)==ceiling(i/500) ){
            cat( ">", date(), " : processed ", i, " of ", nrow( SequenceTable ), " lines.\n" )
        }
    }
    close( outcon )
    ## faster
#    Lines <- apply( SequenceTable, MARGIN=1, function( x ){
#        paste0( ">probe:", chip_type, ":", SequenceTable[ i, 1 ], ";", SequenceTable[ i, 2 ], ":", SequenceTable[ i, 3 ], "; original file:", pgf_file, ";create date ", create_date, "\n", SequenceTable[ i, 4 ], "\n" )
#    })
#    writeLines( Lines, con=outfile, sep="" )
    cat( "finished\n" )

}

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
  fastafile <- paste( chipname, ".fa", sep="" )
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
