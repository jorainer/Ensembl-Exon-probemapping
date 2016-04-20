#######
##
my $script_version="0.3.0";
## generates CDF files from the gene-to-transcript mappings in the database.
## the CDF will contain one probeset per gene containing all probes that target any transcript of the gene. This CDF file can then be used to try to identify alternative splicing (e.g. using FIRMA).
## BE AWARE! Probe sets are defined for each gene defined in Ensembl. A gene that is encoded on a "normal" chromosome and on a patched chromosome will get in Ensembl two different gene ids, thus, we define in such a case two probe sets for the same "gene", which share all or most of their probes. In the analysis the additional gene probe sets can however be excluded by restricting the analysis to chromosomes 1:22, X, Y. This guarantees unique probe sets in the analysis.
## generate the following files:
# 0.3.0: 2013-12-23: this bases on "makeGeneCdfFromDB.pl" but uses the CustomCDF package and the "new" alignment pipeline of probes against cdna and ncrna.
#                  : Important update: probes matching more than one gene (but having only a single genomic alignment) are now included in the CDF. Otherwise all probe sets for pre-miRNAs (defined as genes) and their respective host genes would lack a large amount of their probes.


#######
use IO::File;
use DBI;
use Getopt::Std;
use strict;
use warnings;
use Config::Simple;
use CustomCDF::DBAdaptor;
use CustomCDF::GeneAdaptor;
use CustomCDF::ProbeAdaptor;

##database settings
my $dbname="probemapping_";
my $exon_regions_table="ensembl_exon_regions";
my $probe_to_exon_regions="ensembl_probe_to_exon_regions";
my $username="anonuser";
my $password="";
my $host="localhost";

my $ensembl_version;

## settings for the probe mappings...
my $nr_chrom_hit=1;
my $nr_mm=0;
my $bgmatch=0;		# how often a potential background probe can match (perfect alignment) to the genome.
my $bgmatch_gene=0;	# how often a potential background probe can match (perfect alignment) to a gene.
my $bgmatch_gene_mm1=0;	# how often a potential background probe can match (allowing 1 missmatch) to a gene.
my $bgmatch_mm1=0;	# how often a potential background probe can match (allowing 1 missmatch) to the genome.


## settings for CDF file:
my $GCOS="GC3.0";
my $LE="\n";
my $direction=1;		## are probes sense (1) or antisense (2) to target?
my $NROW=0;

my $cdffile;
my $maxgc=18;
my $min_nr_probes=3;

my $dontstop=0;

my $gff_table="gff_features";

my $define_intronic_probesets=0;


my %option=();
#getopts("i:o:f:",\%option);
getopts("d:e:k:g:m:r:s:Efhi",\%option);
if( $option{ h } ){
    print( "\nmake_gene_cdf_from_db.pl version ".$script_version."\n" );
    print( "Creates a gene-level custom CDF for Affymetrix type microarrays. The CDF will contain one probe set per gene containing all probes that target any transcript/exon of the gene. Note, since in Ensembl the same gene, if encoded also on a patched chromosome, has two distinct gene_ids, a probe set for both is defined in the CDF (with the same probe content). Restriction of the analysis to chromosome 1:22, X and Y avoids the problem of non-unique probe sets.\n\n" );
    print( "usage: make_gene_cdf_from_db.pl -ed[Efghkmrs]\n" );
    print( "parameters:\n" );
    print( "-d database that should be queried. defaults to 'probemapping_'.\n" );
    print( "-e ensembl version (e.g. 49), required.\n" );
    print( "-E if a exon-region file should be created, requires a database table specifying these exon regions (created by the define_nonoverlapping_exon_regions.pl script).\n" );
    print( "-f force generation of CDF and do not break when the probe id does not match the calculated index in the CEL file (just print a warning message).\n" );
    print( "-g maximal number of allowed G and Cs in the probe sequence. defaults to 18 (this will exclude 97183 probes on an Exon array)." );
    print( "-h print this help.\n" );
    print( "-i if set, for each gene, in addition to the exonic probe set, an intronic probe set will be defined, containing all probes matching perfectly to the introns of the gene. An 'in' will be appended to the probe set's name." );
    print( "-k A configuration file with all required settings (database connection settings.\n");
    print( "-m minimal number of probes per probeset. defaults to 3, which means that all probesets with fewer then 3 probes are not included/written to the cdf file.\n" );
    print( "-r nr of rows on the array (has to be specified for the Gene array: 1050).\n" );
    print( "-s (optional) specify the strandnes of the microarray, either sense or antisense for newer generation st microarrays or older 3prime arrays, respectively. If not specified the script will try to guess the strandedness from the name of the microarray platform.\n" );
}
if( $option{ f } ){
    $dontstop=1;
}
if( $option{ i } ){
    $define_intronic_probesets=1;
}
if(  defined( $option{ d } ) ){
	$dbname = $option{ d };
}
if(  defined( $option{g} ) ){
	$maxgc = $option{g};
}
if( defined( $option{m} ) ){
	$min_nr_probes = $option{m};
}
if( defined( $option{r} ) ){
	$NROW = $option{r};
}
if( defined( $option{ s } ) ){
    if( $option{ s } ne "sense" ){
	if( $option{ s } ne "antisense" ){
	    die "-s has to be either sense or antisense!";
	}
    }
    if( $option{ s } eq "antisense" ){
	$direction=2;
    }
}

if( !defined( $option{e} ) ){
  die "Ensembl version has to be submitted with the -e option!\n";
}
else{
  $ensembl_version=$option{e};
}

## the hashes containing the information required to write the CDF
my %GENE2PROBES;	## maps genes to probes: gene_id -> string of probes (; separated)
my %GENE2EXONS;	## maps gene to its exons: gene_id -> string of unique exons (; separated)         ?? DO I NEED THIS ??
my %PROBE2X;	## maps probes to x and y coordinates: probe_id -> x
my %PROBE2Y;	## maps probes to x and y coordinates: probe_id -> y
my %PROBE2SEQUENCE;	## probes to sequence: probe_id -> sequence
my %PROBE2START;	## chromosomel start position of a probe (well, the end is start +25): probe_id -> start
my %PROBE2CHROMOSOME;	## probe_id -> chromosome
my %PROBE2STRAND;	## probe_id -> strand
my %PROBES2GENE;	## mapping probes (; separated) to genes.
my %PROBE2GENES;	## maps single probe id to corresponding gene(s)
my %PROBESINCDF;        ## specifies which probes are actually written to the CDF file. this is important for the
                        ## exon-region file, in order to generate consistent CDF and exon-region files.

#connect to db
my $db_adaptor = CustomCDF::DBAdaptor->new( host=>$host, username=>$username, password=>$password, dbname=>$dbname );
my $dbh = $db_adaptor->dbcon;
my $probe_adaptor=$db_adaptor->get_probe_adaptor();
my $gene_adaptor=$db_adaptor->get_gene_adaptor();

## have to define the platform/chiptype!
my $dummyquery = $dbh->prepare("select platform from probe_information limit 1") or die $dbh->errstr();
$dummyquery->execute() or die $dbh->errstr;
my @dummyresult=$dummyquery->fetchrow_array;
my $platform=$dummyresult[ 0 ];
my $platformlower="\L$platform";
$platformlower=~ s/-//g;
$platformlower=~ s/_//g;

## check whether we can guess the type from the platform.
if( !defined( $option{ s } ) ){
    if( index( $platformlower, "st" )==-1 ){
	$direction=2;
    }
}

my $version_nopoints = $script_version;
$version_nopoints =~ s/\.//g;
$dummyquery->finish;
###### defining the CDF file name.
$cdffile=$platformlower."_".$ensembl_version."_".$version_nopoints."genes";
if( $define_intronic_probesets==1 ){
  $cdffile=$cdffile."inex";
}
$cdffile=$cdffile.".cdf";

my $sth;
my $prep_query;
$sth = $dbh->prepare("select distinct gene_id from probe_to_exon_transcript") or die $dbh->errstr();
$sth->execute() or die $dbh->errstr;

my $infostring = "\nmake_gene_cdf_from_db.pl version $script_version.\n";
$infostring=$infostring."[".localtime()."]\n";
$infostring=$infostring."Ensembl version: $ensembl_version, database name: $dbname.\nmax no. of allowed perfect alignments per probe to the genome: $nr_chrom_hit.\nmax no. of alignments with 1 or 2 missmatches per probe to the genome: $nr_mm.\nmax GC content: $maxgc.\nBackground probes (bg_probes_) definition: no. of genomic alignments: $bgmatch, alignments with 1 mismatch: $bgmatch_mm1, alignments to genes: $bgmatch_gene, alignments with 1 mismatch to genes $bgmatch_gene_mm1.\nminimum number of probes per probeset: $min_nr_probes.\nmicrorray is a ";
if( $direction==1 ){
    $infostring=$infostring."sense target ";
}else{
    $infostring=$infostring."antisense target";
}
$infostring=$infostring."array.\n\ngenerating the following files:\n$cdffile: CDF file that can be converted to binary format using the affxparser Bioconductor package.\n$cdffile.probes.txt: probe sequences, used to generate the probes package for GCRMA preprocessing.\n$cdffile.annotation.txt: the annotation for the probesets defined in the CDF.\n$cdffile.probe-annotation.txt: some additional informations for the probes (index, x, y, chromosome start, strand, targeted genes, exons).\nNote: The CDF will contain one probe set per gene containing all probes that target any transcript/exon of the gene. Note, since in Ensembl the same gene, if encoded also on a patched chromosome, has two distinct gene_ids, a probe set for both is defined in the CDF (with the same probe content). Restriction of the analysis to chromosome 1:22, X and Y avoids the problem of non-unique probe sets.\n";
## skipping: $cdffile.exon.txt: file listing all probes per exon.\n$cdffile.exon.probesets.txt: file defining exon-probesets, combining exons that are targeted by the same probes.\n
if( defined( $option{ r } ) ){
  $infostring=$infostring."nr of rows and columns defined manually: ".$option{r}.".\n";
}

if( $option{ E } ){
  $infostring = $infostring."writing non-overlapping exon region file: ".$platformlower."_".$ensembl_version."_".$version_nopoints."-exon-region-annotation.\n";
}

print $infostring;
open(INFOUT, "> $cdffile.settings");
print INFOUT $infostring;
close(INFOUT);

## calling the function that should do something...
fullVersion_cdna();

print "----------------------------------- fine --------------------------------------\n";


###############################################################################################
# makes a unique array preserving the original ordering
sub makeUnique{
	my @tmp = @_;
	my %seen = ();
	my @uniq = ();
	foreach my $item ( @tmp ){
		push ( @uniq, $item ) unless $seen{ $item }++;
	}
	return @uniq;
}


###############################################################################################
#### should be the full, long version of the function...
sub fullVersion_cdna{
  ## query to fetch intronic probes...
  my $prep_query_intronic=$dbh -> prepare_cached( "select query_pa.probe_id, sequence, probe_information.x, probe_information.y, query_pa.chromosome_start, query_pa.chromosome_strand from (select * from probe_alignments where seq_type='chromosome' and seq_name=? and seq_strand=? and missmatches=0 and start>=? and end<=?) as query_pa join probe_information on (query_pa.probe_id=probe_information.probe_id) where nr_chrom_map<=$nr_chrom_hit and nr_chrom_map_mmall<=$nr_mm and gc_count<=$maxgc and nr_gene_map=0 order by chromosome_start;" );
  my $counter=0;
  ## loop through all genes
  while(my(@row) = $sth->fetchrow_array) {
    my $gene_id = $row[0];
    ## I need:
    ## 1) all probes targeting this gene (inclusively their genomic alignments).
    my @probes = $probe_adaptor->fetch_probes_for_gene( gene_id=>$gene_id, load_genomic_alignments=>1 );
    ### some dummy variables
    my @algs;
    my $alg;
    my @starts;
    my @ends;
    my $probe_id;
    my $dummystring;
    my @genearray = ();
    ## let's loop through the probes.
    foreach my $probe (@probes){
      ## get the alignment of the probe; by their definition we do only have a single genomic alignment for this probes, and no splice junction probes.
      $probe_id = $probe->id();
      @algs = $probe->alignments();
      $alg = $algs[0];
      @starts = $alg->start();
      @ends = $alg->end();
      push( @genearray, $probe_id );  ## this holds then the probe ids for the gene.
      $PROBE2X{ $probe_id } = $probe->index_x;	# do not care if it is overwritten...
      $PROBE2Y{ $probe_id } = $probe->index_y;
      $PROBE2SEQUENCE{ $probe_id } = $probe->sequence;
      $PROBE2START{ $probe_id } = $starts[ 0 ];
      $PROBE2CHROMOSOME{ $probe_id } = $alg->seq_name;
      $PROBE2STRAND{ $probe_id } = $ends[ 0 ];
      if( exists( $PROBE2GENES{ $probe_id } ) ){
	## have also to make shure that i have not duplicated entries...
	$dummystring=$PROBE2GENES{ $probe_id }.";".$gene_id;
	$PROBE2GENES{ $probe_id } = join( ";", makeUnique( split( /;/, $dummystring) ) );
      }
      else{
	$PROBE2GENES{ $probe_id } = $gene_id;
      }
    }
    ## make shure we add something only if there ore some probes!
    if( scalar( @genearray ) > 0 ){
      $GENE2PROBES{ $gene_id } = join( ";", makeUnique( @genearray ) );
    }
    ######
    ## define and add intronic probe set for the current gene.
    ##
    if( $define_intronic_probesets==1 ){
      @genearray = ();
      my $gene_startend_query=$dbh->prepare( "select min(transcript_chrom_start), max(transcript_chrom_end), transcript_chrom_strand, chromosome_name from exon_transcript where gene_id=\'$gene_id\';" ) or die $dbh->errstr();
      $gene_startend_query->execute();
      my @gene_startend_query_result=$gene_startend_query->fetchrow_array;
      $gene_id=$gene_id."intron";
      my $query_strand;
      if( $direction==1 ){
	## sense target:
	if( $gene_startend_query_result[ 2 ]==1 ){
	  $query_strand="-";
	}else{
	  $query_strand="+";
	}
      }
      if( $direction==2 ){
	## sense target:
	if( $gene_startend_query_result[ 2 ]==1 ){
	  $query_strand="+";
	}else{
	  $query_strand="-";
	}
      }
      ## let's make the query...
      $prep_query_intronic -> execute( $gene_startend_query_result[ 3], $query_strand, $gene_startend_query_result[ 0 ], $gene_startend_query_result[ 1 ] ) or die "can not execute query ".$prep_query_intronic -> errstr;
      while(my (@proberow) = $prep_query_intronic->fetchrow_array){
	my $dummystring;
	my $probe_id = $proberow[ 0 ];
	push( @genearray, $probe_id );
	$PROBE2X{ $probe_id } = $proberow[ 2 ];	# do not care if it is overwritten...
	$PROBE2Y{ $probe_id } = $proberow[ 3 ];
	$PROBE2SEQUENCE{ $probe_id } = $proberow[ 1 ];
	$PROBE2START{ $probe_id } = $proberow[ 4 ];
	$PROBE2CHROMOSOME{ $probe_id } = $gene_startend_query_result[ 3 ];
	$PROBE2STRAND{ $probe_id } = $proberow[ 5 ];
	#0.1.6: make a mapping probe to gene(s)
	if( exists( $PROBE2GENES{ $probe_id } ) ){
	  ## have also to make shure that i have not duplicated entries...
	  $dummystring=$PROBE2GENES{ $probe_id }.";".$gene_id;
	  $PROBE2GENES{ $probe_id } = join( ";", makeUnique( split( /;/, $dummystring) ) );
	}
	else{
	  $PROBE2GENES{ $probe_id } = $gene_id;
	}
      }
      ## make shure we add something only if there ore some probes!
      if( scalar( @genearray ) > 0 ){
	$GENE2PROBES{ $gene_id } = join( ";", makeUnique( @genearray ) );
      }
    }
    ## finished with intronic probe set stuff
    $counter++;
    if( ($counter % 1000) == 0 ){
      print localtime()." > processed $counter genes\n";
    }
  }
  $prep_query_intronic->finish;
  print localtime()." > finished querying the database\n";
  addBackgroundProbesOrderGC();
  removeProbesetsBasedOnNrProbes();
  writeTheCDF();
  writeProbeAnnotation();
  if( $option{E} ){
    createProbeRegionFile();
  }
}



###############################################################################################
## add bg probes and collect them to "probesets" based on their GC content.
sub addBackgroundProbesOrderGC{
  print localtime()." > adding potential background probes...";
  my $bgprobequery = $dbh->prepare("select probe_id, x, y, sequence, gc_count from probe_information where nr_chrom_map<=$bgmatch and nr_chrom_map_mm1<=$bgmatch_mm1 and nr_gene_map<=$bgmatch_gene;") or die $dbh->errstr;
  $bgprobequery->execute() or die $dbh->errstr;
  while(my(@proberow) = $bgprobequery->fetchrow_array){
    my $bg_gene_id="bg_probes_gc".$proberow[ 4 ];
    if( exists $GENE2PROBES{ $bg_gene_id } ){
      ## add the probe to the existing entry
      $GENE2PROBES{ $bg_gene_id } = $GENE2PROBES{ $bg_gene_id }.";".$proberow[ 0 ];
    }
    else{
      $GENE2PROBES{ $bg_gene_id } = $proberow[ 0 ];
    }
    $PROBE2GENES{ $proberow[ 0 ] } = $bg_gene_id;
    $PROBE2X{ $proberow[ 0 ] } = $proberow[ 1 ];
    $PROBE2Y{ $proberow[ 0 ] } = $proberow[ 2 ];
    $PROBE2SEQUENCE{ $proberow[ 0 ] } = $proberow[ 3 ];
    $PROBE2START{ $proberow[ 0 ] } = 0;
    $PROBE2CHROMOSOME{ $proberow[ 0 ] } = 0;
    $PROBE2STRAND{ $proberow[ 0 ] } = 0;
  }
  print "finished\n";
}


###############################################################################################
## remove probesets that have fewer than $min_nr_probes probes from all the hashes.
sub removeProbesetsBasedOnNrProbes{
  print localtime()." > removing probesets with less then $min_nr_probes probes\n";
  my $count_removed=0;
  ## saving removed probe sets to a file:
  my $removed_file = $platformlower."_".$ensembl_version."_".$version_nopoints."genes-removed-probesets.txt";
  open(REMOUT, "> $removed_file");
  print REMOUT "gene_id\tno_probes\n";
  foreach my $key (sort keys %GENE2PROBES){
    if( exists( $GENE2PROBES{ $key } ) ){
      my @thedummy = split /;/, $GENE2PROBES{ $key };
      my $nr_probes_ = scalar( @thedummy );
      if( $nr_probes_ < $min_nr_probes ){
	print REMOUT $key."\t".$nr_probes_."\n";
	delete ($GENE2PROBES{ $key });
	#print "deleting $key\n";
	$count_removed++;
      }
    }
    else{
      print "what's with $key?\n";
    }
  }
  close(REMOUT);
  print localtime()." > removed $count_removed probesets\n";
}

###############################################################################################
## try to get the number of rows of the chip based on the biggest y coordinate of a probe. since counting of the rows starts at 0, the maximal y coordinate is incremented by one.
sub getNrOfRows{
  my $rowquery = $dbh->prepare("select max(y) from probe_information") or die $dbh->errstr;
  $rowquery->execute() or die $dbh->errstr;
  my @result=$rowquery->fetchrow_array;
  return ($result[ 0 ] +1);
}

###############################################################################################
## same as getNrOfRows, but using the x coordinates.
sub getNrOfCols{
  my $rowquery = $dbh->prepare("select max(x) from probe_information") or die $dbh->errstr;
  $rowquery->execute() or die $dbh->errstr;
  my @result=$rowquery->fetchrow_array;
  return ($result[ 0 ] +1);
}

###############################################################################################
## calculates the index of a probe in the CEL file. expect to get x and y.
sub getIndex{
  my @xy = @_;
  return ( $xy[ 0 ] + $xy[ 1 ] * $NROW);
}

###############################################################################################
sub complement{
  my $nt = $_[0];
  if( $nt eq "G" ){
    return "C";
  }
  if( $nt eq "C" ){
    return "G";
  }
  if( $nt eq "T" ){
    return "A";
  }
  if( $nt eq "A" ){
    return "T";
  }
}

###############################################################################################
## write the CDF using the hashes... this function is somwhat different fro the function in makeCdfFromDB.pl, since the hashes are different...
## writing additionally also the annotation file, the probe file and a probe information file
sub writeTheCDF{
  my @current_probes;
  my @current_x;
  my @current_y;
  my @current_sequences;
  my $current_probeset_id;
  my $nr_unit=0;
  print localtime()." > writing the CDF file ($cdffile).\n";
  open(OUT, "> $cdffile") or die "can't open output file $cdffile!";
  ## writing the header and chip definition.
  writeCDF();
  writeChip();
  ## writing the annotation...
  open(ANNOTATION, "> $cdffile.annotation.txt");
  print ANNOTATION "probeset_id\tgene_id\tprobe_count\tprobes\n";
  ## writing the probe file...
  open(PROBEFILE, "> $cdffile.probes.txt");
  print PROBEFILE "Probe Set Name\tProbe X\tProbe Y\tProbe Interrogation Position\tProbe Sequence\tTarget Strandedness\n";
  ## ok, now write the Units.
  foreach my $key (sort keys %GENE2PROBES){
    $nr_unit++;
    ## the current probeset id is not so easy to define...
    $current_probeset_id = $ensembl_version.$key;
    @current_probes = split /;/, $GENE2PROBES{ $key };
    my $current_nr_atom=scalar(@current_probes);
    writeUnit( $nr_unit, $current_nr_atom, $current_nr_atom, $nr_unit, 1);
    my $indexstring="";
    ##ok, next write the block.
    print OUT "[Unit".$nr_unit."_Block1]".$LE."Name=$current_probeset_id".$LE."BlockNumber=1".$LE."NumAtoms=".$current_nr_atom.$LE."NumCells=".$current_nr_atom.$LE."StartPosition=0".$LE."StopPosition=".($current_nr_atom-1).$LE."CellHeader=X\tY\tPROBE\tFEAT\tQUAL\tEXPOS\tPOS\tCBASE\tPBASE\tTBASE\tATOM\tINDEX\tCODONIND\tCODON\tREGIONTYPE\tREGION".$LE;
    my $cellcounter = 0;
    my $currentindex=0;
    ## looping through all the probes of the gene and adding it to the CDF file...
    for(my $probe_counter=0;$probe_counter<scalar(@current_probes);$probe_counter++){
      my @nucleotides = split //, $PROBE2SEQUENCE{ $current_probes[ $probe_counter ] };
      my $current_x = $PROBE2X{ $current_probes[ $probe_counter ] };
      my $current_y = $PROBE2Y{ $current_probes[ $probe_counter ] };
      $currentindex = getIndex( $current_x , $current_y );
      if( ($currentindex+1)!=int( $current_probes[ $probe_counter ] ) ){
	if( $dontstop==1 ){
	  print( "WARNING: Estimated number of rows on array seems to be wrong! nrow=".$NROW.", probe id:".$current_probes[ $probe_counter ].", calculated index: ".$currentindex."\n" );
	}
	else{
	  die( "Estimated number of rows on the array seems to be wrong!\nnrow=".$NROW.", probe id:".$current_probes[ $probe_counter ].", calculated index: ".$currentindex."\nYou should manually set the number of rows using the -r option.\n" );
	}
      }
      ## adding the probe to the PROBESINCDF: have to do this, since a probe could be part of several gene-probesets
      ## and therefore probes should not be removed from e.g. PROBE2X hash in the removeProbesetsBasedOnNrProbes
      ##function
      $PROBESINCDF{ $current_probes[ $probe_counter ] } = $current_probes[ $probe_counter ];
      print OUT "Cell".( $probe_counter+1 )."=".$current_x."\t".$current_y."\t"."N"."\tcontrol\t".$current_probeset_id."\t".$probe_counter."\t13\t".complement( $nucleotides[ 12 ] )."\t".$nucleotides[ 12 ]."\t".complement( $nucleotides[ 12 ] )."\t".$probe_counter."\t".$currentindex."\t-1\t-1\t99\t ".$LE;
      ## generating a "index-string" for the annotation table. might be usefull later...
      if( $indexstring ne "" ){
	$indexstring=$indexstring.";".($currentindex+1);
      }
      else{
	$indexstring=($currentindex+1);
      }
      ## writing the probe to the probe file.
      #	print PROBEFILE $current_sequences[ $probe_counter ]."\t".$current_x[ $probe_counter ]."\t".$current_y[ $probe_counter ]."\t".$current_probeset_id."\t".$nr_unit."\tSense\n";
      print PROBEFILE $current_probeset_id."\t".$current_x."\t".$current_y."\t".$nr_unit."\t".$PROBE2SEQUENCE{ $current_probes[ $probe_counter ] }."\tSense\n";
    }
    print OUT "$LE$LE";
    ## write the annotation for the probeset,
    print ANNOTATION $current_probeset_id."\t".$key."\t".$current_nr_atom."\t".$GENE2PROBES{ $key }."\n";
  }
  close(PROBEFILE);
  close(ANNOTATION);
  close(OUT);
}


###############################################################################################
## writing the unit tag.
## expecting the following input: ( probeset_id, numatoms, numcells, unitnumber, numerblocks )
sub writeUnit{
  my @elements = @_;
  print OUT "[Unit$elements[0]]".$LE."Name=NONE".$LE."Direction=$direction".$LE."NumAtoms=$elements[1]".$LE."NumCells=$elements[2]".$LE."UnitNumber=$elements[3]".$LE."UnitType=3".$LE."NumberBlocks=$elements[4]".$LE.$LE;
}


###############################################################################################
## write the CDF block
sub writeCDF{
  print OUT "[CDF]$LE";
  print OUT "Version=$GCOS$LE$LE";
}

###############################################################################################
## write the Chip block
sub writeChip{
  print OUT "[Chip]$LE";
  print OUT "Name=$platform".$LE;
  if( $NROW==0 ){
    $NROW=getNrOfRows();
  }
  my $NCOL=getNrOfCols();
  if( $NROW > $NCOL ){
    $NCOL=$NROW;
  }
  else{
    $NROW=$NCOL;
  }
  print OUT "Rows=".$NROW.$LE;
  print OUT "Cols=".$NCOL.$LE;
  my $hash_length = keys %GENE2PROBES;
  #	my @K = keys %PROBES2TRANSCRIPT;
  print OUT "NumberOfUnits=$hash_length".$LE;
  print OUT "MaxUnit=$hash_length".$LE;
  #	print OUT "MaxUnit=".$K[ ($hash_length - 1) ].$LE;		## according to the Affymetrix documentation, each unit gets its own number, and MaxUnit is then the unit number of the last unit. since we are here using sequential numbering for unit numbers, MaxUnit is equal to NumberOfUnits.
  print OUT "NumQCUnits=0".$LE;
  print OUT "ChipReference= $LE$LE";
}



######################################
#Write the annotation for the individual probes in the CDF #####
# this code was previously in the writeCdf sub, but was removed there, because for probes with more than one
# assignment to a probe set the probe annotation was written for each probe set to the annotation file.
#
sub writeProbeAnnotation{
  ## writing a probe information file...
  open(PROBEINFOFILE, "> $cdffile.probe-annotation.txt");
  print PROBEINFOFILE "probe_id\tcel_index\tx\ty\tchromosome_start\tchromosome_name\tchromosome_strand\tgenes\n";
  foreach my $probe_id (sort keys %PROBESINCDF){
    my $current_x = $PROBE2X{ $probe_id };
    my $current_y = $PROBE2Y{ $probe_id };
    my $currentindex = getIndex( $current_x , $current_y );
    ## well, writing also a probe-information file...
    print PROBEINFOFILE $probe_id."\t".($currentindex+1)."\t".$current_x."\t".$current_y."\t".$PROBE2START{ $probe_id }."\t".$PROBE2CHROMOSOME{ $probe_id }."\t".$PROBE2STRAND{ $probe_id }."\t".$PROBE2GENES{ $probe_id }."\n";
  }
  close( PROBEINFOFILE );
}


######################################
## Write the exon-region file
## this code was taken from define_nonoverlapping_exon_regions.pl
## querying the appropriate database table (exon_probe_regions) for all probes within a exon-regions
## by calling this sub within this script we make shure that the CDF file and the exon-region files are
## consistent, i.e. that the same genes and probes are used/defined/specified in both files.
###
# generate a output file...
sub createProbeRegionFile{
  ## have to define the platform/chiptype!
  my $dummyquery = $dbh->prepare("select platform from probe_information limit 1") or die $dbh->errstr();
  $dummyquery->execute() or die $dbh->errstr;
  my @dummyresult=$dummyquery->fetchrow_array;
  my $platform=$dummyresult[ 0 ];
  my $platformlower="\L$platform";
  $platformlower=~ s/-//g;
  $platformlower=~ s/_//g;
  my $version_nopoints = $script_version;
  $version_nopoints =~ s/\.//g;
  $dummyquery->finish;
  ###### defining the CDF file name.
  my $outfile=$platformlower."_".$ensembl_version."_".$version_nopoints."-exon-region2probes.txt";
  my $outfile_collapsed=$platformlower."_".$ensembl_version."_".$version_nopoints."-exon-region-annotation.txt";
  print localtime()." > writing exon probe region definitions to file $outfile.\n";
  open(OUT, "> $outfile") or die "can't open output file $outfile!";
  open(OUTCOL, "> $outfile_collapsed") or die "can't open output file $outfile_collapsed!";
  print OUT "region_id\tprobe_id\tprobe_chrom_start\n";
  print OUTCOL "region_id\tgene_id\tchromosome_name\tchromosome_strand\texons\tregion_chrom_start\tregion_chrom_end\tprobe_ids\tprobe_count\tpotentially_coding\n";
  ## ok, get all region ids for which we have probes
  my $region_id_query = $dbh->prepare( "select distinct region_id from ".$probe_to_exon_regions ) or die $dbh->errstr;
  $region_id_query->execute() or die $dbh->errstr;
  my $region_counter=0;
  ## prepare a cached query...
  my $prep_query_get_probes_for_region = $dbh->prepare_cached( "select gene_id, chromosome_name, chromosome_strand, exons, probes.region_id, region_chrom_start, region_chrom_end, probe_id, probe_chrom_start, potentially_coding from ( select probe_id, region_id, probe_chrom_start from ".$probe_to_exon_regions." where region_id=(?)) as probes join ".$exon_regions_table." on (probes.region_id=".$exon_regions_table.".region_id) order by probe_chrom_start" ) or die $dbh->errstr;
  while( my @row = $region_id_query->fetchrow_array){
    my $region_id = $row[ 0 ];
    $prep_query_get_probes_for_region -> execute( $region_id ) or die $dbh->errstr;
    my $nr_probes=0;
    my $chromosome_name="";
    my $chromosome_strand;
    my $region_chrom_start;
    my $region_chrom_end;
    my $probe_ids="";
    my $gene_id;
    my $exons;
    my $potentially_coding = 0;
    while( my @regionrow = $prep_query_get_probes_for_region->fetchrow_array ){
      ## jump to next if the probe_id is not in the CDF (means we will not have an expression value for that!)
      if( !exists( $PROBESINCDF{ $regionrow[ 7 ] } ) ){  ## was 6 here... why???
	next;
      }
      if( $nr_probes==0 ){
	$gene_id = $regionrow[ 0 ];
	$chromosome_name = $regionrow[ 1 ];
	$chromosome_strand = $regionrow[ 2 ];
	$exons = $regionrow[ 3 ];
	$region_chrom_start = $regionrow[ 5 ];
	$region_chrom_end = $regionrow[ 6 ];
	$probe_ids = $regionrow[ 7 ];
	$potentially_coding = $regionrow[ 9 ];
	## next exon-region if gene is not present in GENE2PROBES (i.e. in the CDF file)
	if( !exists( $GENE2PROBES{ $gene_id } ) ){
	  last; ## or next...???
	}
      }
      else{
	$probe_ids = $probe_ids.";".$regionrow[ 7 ];
      }
      $nr_probes++;
      # print the region to probe mapping.
      print OUT $region_id."\t".$regionrow[ 7 ]."\t".$regionrow[ 8 ]."\n";
    }
    # print the collapsed exon region file... but only if we have something to write...
    if( $nr_probes > 0 ){
      print OUTCOL $region_id."\t".$gene_id."\t".$chromosome_name."\t".$chromosome_strand."\t".$exons."\t".$region_chrom_start."\t".$region_chrom_end."\t".$probe_ids."\t".$nr_probes."\t".$potentially_coding."\n";
    }
    $region_counter++;
    if( ( $region_counter % 10000 ) == 0 ){
      print localtime()." > processed ".$region_counter." regions\n";
    }
  }
  close( OUT );
  close( OUTCOL );
  $prep_query_get_probes_for_region -> finish;
  print localtime()." > ...finished.\n";
}

################################################################################
##
##    Depracated subs and "old" annotation database with alignment against genome.
##
################################################################################



