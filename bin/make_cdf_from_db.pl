#######
my $script_version="0.3.0";
## generates CDF files from the probe-to-transcript mappings in the database.
# 0.3.0: 2013-12-22: first version bases on makeCDFfromDB.pl but uses the CustomCDF package.
#                  : Important update: probes matching more than one gene (but having only a single genomic alignment) are now included in the CDF. Otherwise all probe sets for pre-miRNAs (defined as genes) and their respective host genes would lack a large amount of their probes.
# 0.3.0: 2014-03-03: fixed a problem with genes encoded on a regular and a patched chromosome:
#                    such a gene has two distinct gene_ids in Ensembl, but can have probe sets
#                    with the same probes. Thus we are now reporting both genes and both chroms
#                    in the gene_id and chromosome_name columns.
#
###
# * TODO: check for HGU133 arrays: probes sets correctly defined based on distance to 3' end?


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

## default database settings
my $dbname="probemapping_";
my $username="anonuser";
my $password="";
my $host="localhost";

my $ensembl_version;

## settings for the probe mappings...
my $nr_chrom_hit=1;
my $nr_mm=0;
my $bgmatch=0;		# how often a potential backgroundprobe can match (perfect alignment) to the genome.
my $bgmatch_gene=0;	# how often a potential backgroundprobe can match (perfect alignment) to a gene.
my $bgmatch_gene_mm1=0;	# how often a potential backgroundprobe can match (allowint 1 mismatch) to a gene.
my $bgmatch_mm1=0;	# how often a potential backgroundprobe can match (allowing 1 missmatch) to the genome.

my $max_distance_3prime=800;
my $check_distance_3prime=0;
my $do_overlapfile=0;

## settings for CDF file:
my $GCOS="GC3.0";
my $LE="\n";
my $direction=1;		## are probes sense (1) or antisense (2) to target?
my $NROW=0;

my $cdffile;
my $maxgc=18;
my $min_nr_probes=3;

my $dontstop=0;

my %option=();
getopts("d:e:g:k:m:r:s:ftph",\%option);

if( $option{ h } ){
## print help... obviously...
    print( "\nmake_cdf_from_db.pl version ".$script_version.".\n" );
    print( "Creates a transcript-level custom CDF for Affymetrix type microarrays. The CDF will contain one probe set for each transcript of a gene containing all probes that target any exon of the transcript. Requires an alignment and annotation database for the respective microarray.\n\n" );
    print( "usage: make_cdf_from_db.pl -d:e:[g:k:m:r:s:ftph]\n" );
    print( "parameters:\n" );
    print( "-d required; database name.\n" );
    print( "-e required; Ensembl version on which the database based.\n" );
    print( "-f force generation of CDF and do not break when the probe id does not match the calculated index in the CEL file (just print a warning message).\n" );
    print( "-g maximal number of allowed G and Cs in the probe sequence. defaults to 18 (this will exclude 97183 probes on an Affymetrix Exon array).\n" );
    print( "-h print this help message.\n" );
    print( "-k A configuration file with all required settings (database connection settings.\n");
    print( "-m minimal number of probes per probeset. defaults to 3, which means that all probesets with fewer then 3 probes are not included/written to the cdf file.\n" );
    print( "-p create a huge tabulator delimited text file (number of probe sets x number of probe sets) that lists the overlap of probe sets in terms of the number of shared probes. The diagonal values contain the numbers of probes per probes.\n" );
    print( "-r nr of rows on the array.\n" );
    print( "-s (optional) specify the strandnes of the microarray, either sense or antisense for newer generation st microarrays or older 3prime arrays, respectively. If not specified the script will try to guess the strandedness from the name of the microarray platform.\n" );
    print( "-t if set probes will only added to a transcripts probe set if their distance to the 3prime end is smaller than $max_distance_3prime.\n\n" );
    exit 0;
}

if( $option{ f } ){
  $dontstop=1;
}
if( $option{ p } ){
  $do_overlapfile=1;
}
if( $option{ t } ){
  $check_distance_3prime=1;
}
if(  defined( $option{g} ) ){
  $maxgc = $option{g};
}
if( defined( $option{m} ) ){
  $min_nr_probes = $option{m};
}
if( defined( $option{ d } ) ){
  $dbname=$option{ d };
}
if( defined( $option{ r } ) ){
  $NROW=$option{ r };
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
  die "Ensembl version has to be submitted with parameter -e\n";
}
else{
  $ensembl_version=$option{e};
}
if( defined( $option{ k } ) ){
  my $cfg = new Config::Simple( $option{ k } );
  $host=$cfg->param( "database.mysql_host" );
  $username=$cfg->param( "database.mysql_user" );
  $password=$cfg->param( "database.mysql_password" );
}

## the hashes containing the information required to write the CDF
my %PROBES2TRANSCRIPT;
my %PROBES2TRANSCRIPTBIOTYPE;
my %PROBES2X;
my %PROBES2Y;
my %PROBES2SEQUENCE;
my %PROBES2GENE;
### huge new hashes....
my %PROBEINFORMATION;
my %PROBEINFORMATION2TRANSCRIPT;
my %PROBEINFORMATION2GENE;

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

my $is_cdna_db = 0;
if( index( $dbname, "cdna" )==-1 ){
  ## well, have an "old" style db.
}else{
  $is_cdna_db = 1;
}

#$dummyquery->finish;
my $version_nopoints = $script_version;
$version_nopoints =~ s/\.//g;
###### defining the CDF file name.
$cdffile=$platformlower."_".$ensembl_version."_".$version_nopoints.".cdf";

my $overlapfile = $platformlower."_".$ensembl_version."_".$version_nopoints."-probeset-overlap.txt";

## workflow:
# 1) get a unque list of transcripts from the database.
# 2) for each transcript get all probes from the database..

# 1) this gets a list of all transcripts that are targeted by at least on probe.
my $sth = $dbh->prepare("select distinct transcript_id from probe_to_exon_transcript") or die $dbh->errstr();
$sth->execute() or die $dbh->errstr;

my $prep_query;
my $prep_query_transcript;
my $prep_query_length_fw;
my $prep_query_length_rv;
my $prep_query_get_exon;

## ok, now we really start...
my $infostring = "\nmake_cdf_from_db.pl version $script_version.\n";
$infostring=$infostring."[".localtime()."]\n";
$infostring=$infostring."Ensembl version: $ensembl_version, database name: $dbname.\nmax no. of allowed perfect alignments per probe to the genome: $nr_chrom_hit.\nmax no. of alignments with 1 missmatch per probe to the genome: $nr_mm (2 or more missmatches are allowed).\nmax GC content: $maxgc.\nBackground probes (bg_probes_) definition: no. of genomic alignments: $bgmatch, alignments with 1 mismatch: $bgmatch_mm1, alignments to genes: $bgmatch_gene, alignments to gene(s) with one mismatch $bgmatch_gene_mm1.\nminimum number of probes per probeset: $min_nr_probes.\nmicroarray is a ";
if( $direction==1 ){
    $infostring=$infostring."sense target array.\n";
}else{
    $infostring=$infostring."antisense target array.\n";
}
if( $check_distance_3prime==1 ){
    $infostring=$infostring."probes sets contain only probes that target the mRNA within ".$max_distance_3prime." of the transcripts 3prime end.\n";
}
$infostring=$infostring."\ngenerating the following files:\n$cdffile: CDF file that can be converted to binary format using the affxparser Bioconductor package.\n$cdffile.probes.txt: probe sequences, used to generate the probes package for GCRMA preprocessing.\n$cdffile.annotation.txt: the annotation for the probesets defined in the CDF.\n\n";
if( $is_cdna_db==0 ){
  $infostring="Note: the database you are querying is based on a deprecated alignment and annotation pipeline; you should consider to use the newer cdna_alignment.pl perl script instead.\n".$infostring;
}
print $infostring;

open(INFOUT, "> $cdffile.settings");
print INFOUT $infostring;
close(INFOUT);

## calling the function that should do something...
if( $is_cdna_db==1 ){
  ##> remove that once I've made sure that it is OK.
  if( $check_distance_3prime==1 ){
    print( "\nWARNING! This option has not been tested extensively!\n" );
  }
  ##<
  ## that's the new cdna and ncrna alignment approach
  fullVersion_cdna();
  undef $probe_adaptor;
  undef $gene_adaptor;
}else{
  print( "\nDEPRECATED!\nPlease query a database that was generated using the cdna_alignment.pl script!\n" );
  ## that's the old genome alignment approach
  ## it is ok to use nr_chrom_map<=1 for examle, since no probes with no alignment are related to trascripts.
  my $prep_query=$dbh -> prepare_cached("select transcript_id, probe_alignments.probe_id, sequence, probe_information.x, probe_information.y, gene_id, exon_id, probe_alignments.chromosome_start, probe_alignments.chromosome_end from (select transcript_id, probes_fk, gene_id, exon_id from probe_to_exon_transcript where transcript_id=?) as queryt join probe_alignments on (probes_fk=probes_pk) join probe_information on (probe_alignments.probe_id=probe_information.probe_id) where nr_chrom_map<=$nr_chrom_hit and nr_chrom_map_mm1<=$nr_mm and gc_count<=$maxgc;");

  my $prep_query_transcript = $dbh -> prepare_cached( "select distinct transcript_biotype, transcript_chrom_strand from exon_transcript where transcript_id=?" ) or die $dbh-> errstr();

  my $prep_query_length_fw = $dbh -> prepare_cached( "select sum(exon_chrom_end-exon_chrom_start) from exon_transcript where transcript_id=? and exon_chrom_start > ?" ) or die $dbh-> errstr();
  my $prep_query_length_rv = $dbh -> prepare_cached( "select sum(exon_chrom_end-exon_chrom_start) from exon_transcript where transcript_id=? and exon_chrom_end < ?" ) or die $dbh-> errstr();
  my $prep_query_get_exon=$dbh->prepare_cached( "select exon_chrom_start,exon_chrom_end from exon_transcript where exon_id=?" ) or die $dbh->errstr();
  fullVersion();
}


##
if( $do_overlapfile==1 ){
  getOverlapOfProbesets();
}

if( $is_cdna_db==0 ){
  ## cleaning up
  $prep_query->finish();
  $prep_query_transcript->finish();
  $prep_query_get_exon->finish();
  $prep_query_length_fw->finish();
  $prep_query_length_rv->finish();
  $dbh->disconnect();
}
print "----------------------------------- fine --------------------------------------\n";

## thats's the "new" version for the alignments against cdna and ncrna.
sub fullVersion_cdna{
  my $counter=0;
  my $add_probe=1;
  while(my(@row) = $sth->fetchrow_array) {
    my $transcript_id = $row[0];
    my $transcript = $gene_adaptor->fetch_transcript( transcript_id=>$transcript_id, load_exons=>1 );
    my $gene = $transcript->fetch_gene;
    my $gene_id = $gene->gene_id;
    my @probes;
    if( $check_distance_3prime==1 ){
      ## 2) now get all probes with their alignments (to cdna) for the transcript.
      ## we just need the alignments for the 3' arrays.
      @probes = $probe_adaptor->fetch_probes_for_transcript( transcript_id=>$transcript_id , load_alignments=>1, seq_type=>"='cdna'");
    }else{
      ## should be faster.
      @probes = $probe_adaptor->fetch_probes_for_transcript( transcript_id=>$transcript_id );
    }
    my $transcript_biotype=$transcript->transcript_biotype;
    my $probestring="";
    my $xstring="";
    my $ystring="";
    my $sequencestring="";
    my $transcript_chrom_strand=$transcript->strand+0;
    my $transcript_length=$transcript->length;
    ## now loop through the probes...
    foreach my $probe (@probes){
      ## do not add the probe if it has more than $maxgc G and C nucleotides.
      $add_probe=1;
      if( $probe->gc_count>$maxgc ){
	$add_probe=0;
      }
      ## for microarrays with random primed targets it is straight forward, just add all probes.
      ## for oligoDT primed targets its complicated since there is a 3' bias and it makes only
      ## sense to add probes to a probe set if they are not too far from the 3' end.
      if( $check_distance_3prime==1 ){
	## have a 3' array, thus checking if distance to 3' end is less than the max.
	my @algs =$ probe->alignments();
	my $alg = $algs[ 0 ]; ## assuming a single alignment, otherwise we do have a problem!
	my @starts = $alg->start;
	$add_probe = 0;
	if( ( $transcript_length-$starts[ 0 ] ) < $max_distance_3prime ){
	  $add_probe=1;
	}
      }
      ## ok and now move on.
      if( $add_probe==1 ){
	if( $probestring ne ""){
	  $probestring=$probestring.";".$probe->id();
	  $xstring = $xstring.";".$probe->index_x();
	  $ystring = $ystring.";".$probe->index_y();
	  $sequencestring = $sequencestring.";".$probe->sequence();
	}
	else{
	  $probestring=$probestring.$probe->id();
	  $xstring = $xstring.$probe->index_x();
	  $ystring = $ystring.$probe->index_y();
	  $sequencestring = $sequencestring.$probe->sequence;
	}
      }
    }
    if( $probestring ne ""){
      ## add all values to the hashes.
      if( exists $PROBES2TRANSCRIPT{ $probestring } ){
	## we do have already a probe set with the same probes.
	my $dummy = $PROBES2TRANSCRIPT{ $probestring }.";".$transcript_id;
	$PROBES2TRANSCRIPT{$probestring} = $dummy;
	$dummy = $PROBES2TRANSCRIPTBIOTYPE{ $probestring }.";".$transcript_biotype;
	$PROBES2TRANSCRIPTBIOTYPE{ $probestring } = $dummy;
	my @genearray = split( /;/, $PROBES2GENE{ $probestring } );
	push ( @genearray, $gene_id );
	$PROBES2GENE{ $probestring } = join( ";", makeUnique( @genearray ) );
      }
      else{
	$PROBES2TRANSCRIPT{$probestring} = $transcript_id;
	$PROBES2X{$probestring} = $xstring;
	$PROBES2Y{$probestring} = $ystring;
	$PROBES2SEQUENCE{$probestring} = $sequencestring;
	$PROBES2GENE{$probestring} = $gene_id;
	$PROBES2TRANSCRIPTBIOTYPE{$probestring} = $transcript_biotype;
      }
    }
    $counter++;
    if( ($counter % 5000) == 0 ){
      print localtime()." > processed $counter transcripts\n";
    }
  }
  print localtime()." > finished querying the database\n";
  addBackgroundProbesOrderGC();
  removeProbesetsBasedOnNrProbes();
  writeTheCDF();
}

###########
## add bg probes and collect them to "probesets" based on their GC content.
sub addBackgroundProbesOrderGC{
  print localtime()." > adding potential background probes...";
  my $bgprobequery;
  if( $is_cdna_db==1 ){
    $bgprobequery = $dbh->prepare("select probe_id, x, y, sequence, gc_count from probe_information where nr_chrom_map<=$bgmatch and nr_chrom_map_mm1<=$bgmatch_mm1 and nr_gene_map<=$bgmatch_gene and nr_gene_map_mm<=$bgmatch_gene_mm1;") or die $dbh->errstr;
  }else{
    $bgprobequery = $dbh->prepare("select probe_id, x, y, sequence, gc_count from probe_information where nr_chrom_map<=$bgmatch and nr_chrom_map_mm1<=$bgmatch_mm1 and nr_gene_map<=$bgmatch_gene;") or die $dbh->errstr;
  }
  $bgprobequery->execute() or die $dbh->errstr;
  while(my(@proberow) = $bgprobequery->fetchrow_array){
    my $bg_probe_id="bg_probes_gc".$proberow[ 4 ];
    if( exists $PROBES2TRANSCRIPT{ $bg_probe_id } ){
      ## add it to the existing entry
      #my $dummy = $PROBES2TRANSCRIPT{ $probestring }.";".$transcript_id;
      #$PROBES2TRANSCRIPT{$probestring} = $dummy;
      $PROBES2X{ $bg_probe_id } = $PROBES2X{ $bg_probe_id }.";".$proberow[ 1 ];
      $PROBES2Y{ $bg_probe_id } = $PROBES2Y{ $bg_probe_id }.";".$proberow[ 2 ];
      $PROBES2SEQUENCE{ $bg_probe_id } = $PROBES2SEQUENCE{ $bg_probe_id }.";".$proberow[ 3 ];
    }
    else{
      $PROBES2TRANSCRIPT{ $bg_probe_id } = $bg_probe_id;
      $PROBES2TRANSCRIPTBIOTYPE{ $bg_probe_id } = "exon_background_probes";
      $PROBES2X{ $bg_probe_id } = $proberow[ 1 ];
      $PROBES2Y{ $bg_probe_id } = $proberow[ 2 ];
      $PROBES2SEQUENCE{ $bg_probe_id } = $proberow[ 3 ];
      $PROBES2GENE{ $bg_probe_id } = "BGPROBE_GC".$proberow[ 4 ];
    }
  }
  print "finished\n";
}

############
## remove probesets that have fewer than $min_nr_probes probes from all the hashes.
sub removeProbesetsBasedOnNrProbes{
  print localtime()." > removing probesets with less then $min_nr_probes probes\n";
  ## saving removed probe sets to a file:
  my $removed_file=$platformlower."_".$ensembl_version."_".$version_nopoints."-removed-probesets.txt";
  open(REMOUT, "> $removed_file");
  print REMOUT "gene_id\ttranscript_id\tno_probes\tno_probesets_kept_for_gene\n";
  my $count_removed=0;
  my %REMPROBES2GENE;
  my %REMPROBES2TRANSCRIPT;
  foreach my $key (sort keys %PROBES2X){
    my @thedummy = split( /;/, $PROBES2X{ $key } );
    my $nr_probes_ = scalar( @thedummy );
    if( $nr_probes_ < $min_nr_probes ){
      $REMPROBES2GENE{ $key } = $PROBES2GENE{ $key };
      $REMPROBES2TRANSCRIPT{ $key } = $PROBES2TRANSCRIPT{ $key };
      delete ($PROBES2X{ $key });
      delete ($PROBES2Y{ $key });
      delete ($PROBES2TRANSCRIPT{ $key });
      delete ($PROBES2TRANSCRIPTBIOTYPE{ $key });
      delete ($PROBES2GENE{ $key });
      delete ($PROBES2SEQUENCE{ $key });
      $count_removed++;
    }
  }
  ## make a mapping gene to no. of probe sets.
  my %KEPTGENES2PROBESETNO=();
  foreach my $key (sort keys %PROBES2X){
    my $the_gene = $PROBES2GENE{ $key };
    if( exists $KEPTGENES2PROBESETNO{ $the_gene } ){
      $KEPTGENES2PROBESETNO{ $the_gene } = $KEPTGENES2PROBESETNO{ $the_gene } + 1;
    }else{
      $KEPTGENES2PROBESETNO{ $the_gene } = 1;
    }
  }
  foreach my $key (sort keys %REMPROBES2GENE ){
    my $probesets_left = 0;
    my $the_gene = $REMPROBES2GENE{ $key };
    if( exists $KEPTGENES2PROBESETNO{ $the_gene } ){
      $probesets_left = $KEPTGENES2PROBESETNO{ $the_gene };
    }
    my @thedummy = split( /;/, $key );
    my $nr_probes_ = scalar( @thedummy );
    print REMOUT $the_gene."\t".$REMPROBES2TRANSCRIPT{ $key }."\t".$nr_probes_."\t".$probesets_left."\n";
  }
  close(REMOUT);
  print localtime()." > removed $count_removed probesets\n";
}

########
## try to get the number of rows of the chip based on the biggest y coordinate of a probe. since counting of the rows starts at 0, the maximal y coordinate is incremented by one.
sub getNrOfRows{
  my $rowquery = $dbh->prepare("select max(y) from probe_information") or die $dbh->errstr;
  $rowquery->execute() or die $dbh->errstr;
  my @result=$rowquery->fetchrow_array;
  return ($result[ 0 ] +1);
}

########
## same as getNrOfRows, but using the x coordinates.
sub getNrOfCols{
  my $rowquery = $dbh->prepare("select max(x) from probe_information") or die $dbh->errstr;
  $rowquery->execute() or die $dbh->errstr;
  my @result=$rowquery->fetchrow_array;
  return ($result[ 0 ] +1);
}

## calculates the index of a probe in the CEL file. expect to get x and y.
sub getIndex{
  my @xy = @_;
  return ( $xy[ 0 ] + $xy[ 1 ] * $NROW);
}

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

##################
## write the CDF using the hashes...
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
  print ANNOTATION "probeset_id\ttranscript_id\tgene_id\tprobe_count\tprobes\ttranscript_biotype\n";
  #	print ANNOTATION "probeset_id\ttranscript_id\tgene_id\tnr_probes\tprobes\tindices\n";

  ## writing the probe file...
  open(PROBEFILE, "> $cdffile.probes.txt");
  print PROBEFILE "Probe Set Name\tProbe X\tProbe Y\tProbe Interrogation Position\tProbe Sequence\tTarget Strandedness\n";

  ## ok, now write the Units.
  foreach my $key (sort keys %PROBES2X){
    $nr_unit++;
    ## the current probeset id is not so easy to define...
    my @alltranscripts = split( /;/, $PROBES2TRANSCRIPT{ $key });
    $current_probeset_id = $ensembl_version.$alltranscripts[0];
    #		$current_probeset_id = $nr_unit."_".$ensembl_version;
    @current_probes = split /;/,$key;
    @current_x = split /;/, $PROBES2X{ $key };
    @current_y = split /;/, $PROBES2Y{ $key };
    @current_sequences = split /;/, $PROBES2SEQUENCE{ $key };
    my $current_nr_atom=scalar(@current_x);
    writeUnit( $nr_unit, $current_nr_atom, $current_nr_atom, $nr_unit, 1);

    my $indexstring="";
    ##ok, next write the block.
    print OUT "[Unit".$nr_unit."_Block1]".$LE."Name=$current_probeset_id".$LE."BlockNumber=1".$LE."NumAtoms=".$current_nr_atom.$LE."NumCells=".$current_nr_atom.$LE."StartPosition=0".$LE."StopPosition=".($current_nr_atom-1).$LE."CellHeader=X\tY\tPROBE\tFEAT\tQUAL\tEXPOS\tPOS\tCBASE\tPBASE\tTBASE\tATOM\tINDEX\tCODONIND\tCODON\tREGIONTYPE\tREGION".$LE;
    my $cellcounter = 0;
    my $currentindex=0;
    ## looping through the x indices and not the probe ids, because the probe ids are not submitted for the background probes in the $key attribute.
    for(my $probe_counter=0;$probe_counter<scalar(@current_x);$probe_counter++){
      my @nucleotides = split //,$current_sequences[ $probe_counter ];
      $currentindex = getIndex( $current_x[ $probe_counter ], $current_y[ $probe_counter ] );
      ## paranoia!!!!
      ## check if the probe id is identical to the cel index calculated
      ## based on x and y
      if( $PROBES2GENE{ $key } !~ /^BGPROBE/ ){
	if( ($currentindex+1)!=int( $current_probes[ $probe_counter ] ) ){
	  if( $dontstop==1 ){
	    print( "WARNING: Estimated number of rows on array seems to be wrong! nrow=".$NROW.", probe id:".$current_probes[ $probe_counter ].", calculated index: ".$currentindex."\n" );
	  }
	  else{
	    die( "Estimated number of rows on the array seems to be wrong!\nnrow=".$NROW.", probe id:".$current_probes[ $probe_counter ].", calculated index: ".$currentindex."\nYou should manually set the number of rows using the -r option.\n" );
	  }
	}
      }
      print OUT "Cell".( $probe_counter+1 )."=".$current_x[ $probe_counter ]."\t".$current_y[ $probe_counter ]."\t"."N"."\tcontrol\t".$current_probeset_id."\t".$probe_counter."\t13\t".complement( $nucleotides[ 12 ] )."\t".$nucleotides[ 12 ]."\t".complement( $nucleotides[ 12 ] )."\t".$probe_counter."\t".$currentindex."\t-1\t-1\t99\t ".$LE;

      ## generating a "index-string" for the annotation table. might be usefull later...
      if( $indexstring ne "" ){
	$indexstring=$indexstring.";".($currentindex+1);
      }
      else{
	$indexstring=($currentindex+1);
      }

      ## writing the probe to the probe file.
      #			print PROBEFILE $current_sequences[ $probe_counter ]."\t".$current_x[ $probe_counter ]."\t".$current_y[ $probe_counter ]."\t".$current_probeset_id."\t".$nr_unit."\tSense\n";
      print PROBEFILE $current_probeset_id."\t".$current_x[ $probe_counter ]."\t".$current_y[ $probe_counter ]."\t".$nr_unit."\t".$current_sequences[ $probe_counter ]."\tSense\n";
    }
    print OUT "$LE$LE";

    ## write the annotation for the probeset,
    print ANNOTATION $current_probeset_id."\t".$PROBES2TRANSCRIPT{ $key }."\t".$PROBES2GENE{ $key }."\t".$current_nr_atom."\t".$key."\t".$PROBES2TRANSCRIPTBIOTYPE{ $key }."\n";

  }

  close(PROBEFILE);
  close(ANNOTATION);
  close(OUT);
}


#####
## writing the unit tag.
## expecting the following input: ( probeset_id, numatoms, numcells, unitnumber, numerblocks )
sub writeUnit{
  my @elements = @_;
  print OUT "[Unit$elements[0]]".$LE."Name=NONE".$LE."Direction=$direction".$LE."NumAtoms=$elements[1]".$LE."NumCells=$elements[2]".$LE."UnitNumber=$elements[3]".$LE."UnitType=3".$LE."NumberBlocks=$elements[4]".$LE.$LE;
}


## write the CDF block
sub writeCDF{
  print OUT "[CDF]$LE";
  print OUT "Version=$GCOS$LE$LE";
}

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
  my $hash_length = keys %PROBES2X;
  #	my @K = keys %PROBES2TRANSCRIPT;
  print OUT "NumberOfUnits=$hash_length".$LE;
  print OUT "MaxUnit=$hash_length".$LE;
  #	print OUT "MaxUnit=".$K[ ($hash_length - 1) ].$LE;		## according to the Affymetrix documentation, each unit gets its own number, and MaxUnit is then the unit number of the last unit. since we are here using sequential numbering for unit numbers, MaxUnit is equal to NumberOfUnits.
  print OUT "NumQCUnits=0".$LE;
  print OUT "ChipReference= $LE$LE";
}

sub getOverlapOfProbesets{
## loop over probe sets, i.e. keys in PROBE2X
## with two idx, i and j
## ... how to make a for loop with an index in perl...
## probe2x is a hash!
  print localtime()." > writing the probe set overlap file ($overlapfile).\n";
  open(OUT, "> $overlapfile") or die "can't open output file $overlapfile!";
  my @keys = sort keys %PROBES2X;
  ## mmmh, unefficient code below...
  my $prefix="";
  my $dummyout="";
  foreach my $current_key (@keys){
    my @alltranscripts = split( /;/, $PROBES2TRANSCRIPT{ $current_key });
    $dummyout=$dummyout.$prefix.$ensembl_version.$alltranscripts[0];
    $prefix="\t";
  }
  print OUT $dummyout."\n";
  ## ok that was the header
  ##
  ## print OUT join( "\t", @keys )."\n";
  my $i; ## rows
  my $j; ## cols
  for ($i=0; $i < scalar (@keys); $i++ ){
    my $keya = $keys[ $i ];
    my @a = split( /;/, $keya );
    my $outstring="";
    if( $i > 0 ){
      $outstring= ( "\t" x $i );
    }
    ## at least that one should get increasingly faster.
    for( $j=$i; $j < scalar (@keys); $j++ ){
      my $keyb = $keys[ $j ];
      my @b = split( /;/, $keyb );
      my %union;
      my %intersect;
      foreach my $e (@a, @b){
	$union{ $e }++ && $intersect{ $e }++;
      }
      if( $j > $i ){
	$outstring=$outstring."\t";
      }
      $outstring=$outstring.scalar( keys %intersect );
    }
    print OUT $outstring."\n";
    if( ($i % 1000) == 0 ){
      print localtime()." > processed $i of in total ".scalar @keys." probe sets.\n";
    }
  }
  close( OUT );
}



################################################################################
##
##    Depracated subs and "old" annotation database with alignment against genome.
##
################################################################################

sub fullVersion{
    my $counter=0;
    my $add_probe=1;
    while(my(@row) = $sth->fetchrow_array) {
	my $transcript_id = $row[0];
	my $gene_id = "";
	## 2) now get all probes targetting this transcript.
	$prep_query -> execute( $transcript_id ) or die "can not execute query ".$prep_query -> errstr;
	$prep_query_transcript -> execute( $transcript_id ) or die "can not execute query ".$prep_query_transcript -> errstr;
	## assuming i get only one entry...
	my (@dummyresult) = $prep_query_transcript->fetchrow_array;
	my $transcript_biotype=$dummyresult[ 0 ];
	my $probestring="";
	my $xstring="";
	my $ystring="";
	my $sequencestring="";
	my $transcript_chrom_strand=$dummyresult[ 1 ]+0;  ## get strand and make shure its a number...
	while(my (@proberow) = $prep_query->fetchrow_array){
	    ## for microarrays with random primed targets it is straight forward, just add all probes.
	    ## for oligoDT primed targets its complicated since there is a 3' bias and it makes only
	    ## sense to add probes to a probe set if they are not too far from the 3' end.
	    if( $check_distance_3prime==1 ){
		my $distance_3prime=getDistance3Prime( $transcript_id, $transcript_chrom_strand, $proberow[ 6 ], $proberow[ 7 ], $proberow[ 8 ] );
		if( $distance_3prime > $max_distance_3prime ){
		    $add_probe=0;
		}else{
		    $add_probe=1;
		}
	    }
	    if( $add_probe==1 ){
		if( $probestring ne ""){
		    $probestring=$probestring.";".$proberow[1];
		    $xstring = $xstring.";".$proberow[ 3 ];
		    $ystring = $ystring.";".$proberow[ 4 ];
		    $sequencestring = $sequencestring.";".$proberow[ 2 ];
		}
		else{
		    $probestring=$probestring.$proberow[1];
		    $xstring = $xstring.$proberow[ 3 ];
		    $ystring = $ystring.$proberow[ 4 ];
		    $sequencestring = $sequencestring.$proberow[ 2 ];
		    $gene_id = $proberow[ 5 ];
		}
	    }
	}
	if( $probestring ne ""){
	    ## add all values to the hashes.
	    if( exists $PROBES2TRANSCRIPT{ $probestring } ){
		my $dummy = $PROBES2TRANSCRIPT{ $probestring }.";".$transcript_id;
		$PROBES2TRANSCRIPT{$probestring} = $dummy;
		$dummy = $PROBES2TRANSCRIPTBIOTYPE{ $probestring }.";".$transcript_biotype;
		$PROBES2TRANSCRIPTBIOTYPE{ $probestring } = $dummy;
	    }
	    else{
		$PROBES2TRANSCRIPT{$probestring} = $transcript_id;
		$PROBES2X{$probestring} = $xstring;
		$PROBES2Y{$probestring} = $ystring;
		$PROBES2SEQUENCE{$probestring} = $sequencestring;
		$PROBES2GENE{$probestring} = $gene_id;
		$PROBES2TRANSCRIPTBIOTYPE{$probestring} = $transcript_biotype;
	    }
	}
	$counter++;
	if( ($counter % 5000) == 0 ){
	    print localtime()." > processed $counter transcripts\n";
	}
    }
    print "finished querying the database\n";
    if( $option{a} ){
	addProbesetsFromGFF();	## put this after removeProbesetsBasedOnNrProbes if the no probesets should be removed for this set...
    }
    addBackgroundProbesOrderGC();
    removeProbesetsBasedOnNrProbes();
    writeTheCDF();
}

#close(PROBEOUT);

###########
## adding potential "background probes", i.e. probes not matching to the genome.
sub addBackgroundProbes{
	print "adding potential background probes...";
	my $bgprobequery = $dbh->prepare("select probe_id, x, y, sequence from probe_information where nr_chrom_map<=$bgmatch and nr_chrom_map_mm1<=$bgmatch_mm1 and nr_gene_map<=$bgmatch_gene;") or die $dbh->errstr;
	$bgprobequery->execute() or die $dbh->errstr;

	while(my(@proberow) = $bgprobequery->fetchrow_array){
		$PROBES2TRANSCRIPT{$proberow[ 0 ]} = "bg_".$proberow[ 0 ];
		$PROBES2TRANSCRIPTBIOTYPE{$proberow[ 0 ]} = "exon_background_probes";
		$PROBES2X{$proberow[ 0 ]} = $proberow[ 1 ];
		$PROBES2Y{$proberow[ 0 ]} = $proberow[ 2 ];
		$PROBES2SEQUENCE{$proberow[ 0 ]} = $proberow[ 3 ];
		$PROBES2GENE{ $proberow[ 0 ] } = "BGPROBE";
	}
	print "finished\n";
}



##############################
## define additional probesets from the information stored in the gff table.
sub addProbesetsFromGFF{
  my $gff_table;
	print "adding additional probesets form the $gff_table table...";
	my $prep_query2=$dbh -> prepare_cached("select probeset_id, probes.probe_id, sequence, probe_information.x, probe_information.y, name, transcript from (select transcript, probes_fk, probeset_id, name from $gff_table where probeset_id=?) as queryt join probes on (probes_fk=probes_pk) join probe_information on (probes.probe_id=probe_information.probe_id) where nr_chrom_map<=$nr_chrom_hit and nr_chrom_map_mm1<=$nr_mm and gc_count<=$maxgc;");

	my $gffquery = $dbh->prepare("select distinct probeset_id from $gff_table") or die $dbh->errstr();
	$gffquery->execute() or die $dbh->errstr;

	while(my(@row) = $gffquery->fetchrow_array) {
		my $transcript_id = $row[0];
		my $gene_id = "";
		## 2) now get all probes targetting this transcript.
		$prep_query2 -> execute( $transcript_id ) or die "can not execute query ".$prep_query -> errstr;
		my $probestring="";
		my $xstring="";
		my $ystring="";
		my $sequencestring="";
		while(my (@proberow) = $prep_query2->fetchrow_array){
			## would be faster and less memory cconsuming if i would just write the results line by line...
			if( $probestring ne ""){
				$probestring=$probestring.";".$proberow[1];
				$xstring = $xstring.";".$proberow[ 3 ];
				$ystring = $ystring.";".$proberow[ 4 ];
				$sequencestring = $sequencestring.";".$proberow[ 2 ];
			}
			else{
				$probestring=$probestring.$proberow[1];
				$xstring = $xstring.$proberow[ 3 ];
				$ystring = $ystring.$proberow[ 4 ];
				$sequencestring = $sequencestring.$proberow[ 2 ];
				$gene_id = $proberow[ 5 ];
			}

#			### new 25.07.2008: store probe information to a file. may be used later to detect splice variants...
#			print PROBEOUT $proberow[ 1 ]."\t".(getIndex( $proberow[ 3 ], $proberow[ 4 ])+1)."\t".$proberow[ 3 ]."\t".$proberow[ 4 ]."\t".$proberow[ 9 ]."\t".$proberow[ 10 ]."\t".$proberow[ 5 ]."\t".$proberow[ 0 ]."\t".$proberow[ 6 ]."\t".$proberow[ 7 ]."\t".$proberow[ 8 ]."\t".$proberow[ 11 ]."\n";

		}
		if( $probestring ne ""){
			## add all values to the hashes.
			if( exists $PROBES2TRANSCRIPT{ $probestring } ){
				my $dummy = $PROBES2TRANSCRIPT{ $probestring }.";".$transcript_id;
				$PROBES2TRANSCRIPT{$probestring} = $dummy;
				$dummy = $PROBES2TRANSCRIPTBIOTYPE{ $probestring }.";".$gff_table;
				$PROBES2TRANSCRIPTBIOTYPE{ $probestring } = $dummy;
			}
			else{
				$PROBES2TRANSCRIPT{$probestring} = $transcript_id;
				$PROBES2TRANSCRIPTBIOTYPE{$probestring} = $gff_table;
				$PROBES2X{$probestring} = $xstring;
				$PROBES2Y{$probestring} = $ystring;
				$PROBES2SEQUENCE{$probestring} = $sequencestring;
				$PROBES2GENE{$probestring} = $gene_id;
			}
		}
	}

	print "...finished.\n";
	$prep_query2->finish();
}



## this function does simply write all information to the CDF after each fetched row from the database. there is no checking if some transcripts are targeted with exactly the same probes and merging them.
## DOES NOT WORK!!! for the CDF i have to know beforehand how many probesets i will have!
sub quickCDF{
	open(OUT, "> $option{o}") or die "can't open output file $option{o}!";

	while(my(@row) = $sth->fetchrow_array) {
		my $transcript_id = $row[0];
		## 2) now get all probes targetting this transcript.
		$prep_query -> execute( $transcript_id ) or die "can not execute query ".$prep_query -> errstr;
		my $probestring="";
		my $xstring="";
		my $ystring="";
		my $sequencestring="";
		while(my (@proberow) = $prep_query->fetchrow_array){
			## would be faster and less memory cconsuming if i would just write the results line by line...
			if( $prep_query ne ""){
				$probestring=$probestring.";".$proberow[1];
				$xstring = $xstring.";".$proberow[ 3 ];
				$ystring = $ystring.";".$proberow[ 4 ];
				$sequencestring = $sequencestring.";".$proberow[ 2 ];
			}
			else{
				$probestring=$probestring.$proberow[1];
				$xstring = $xstring.";".$proberow[ 3 ];
				$ystring = $ystring.";".$proberow[ 4 ];
				$sequencestring = $sequencestring.$proberow[ 2 ];
			}
		}
		if( $probestring ne ""){
			## add all values to the hashes.
			if( exists $PROBES2TRANSCRIPT{ $probestring } ){
				my $dummy = $PROBES2TRANSCRIPT{ $probestring }.";".$transcript_id;
				$PROBES2TRANSCRIPT{$probestring} = $dummy;
			}
			else{
				$PROBES2TRANSCRIPT{$probestring} = $transcript_id;
				$PROBES2X{$probestring} = $xstring;
				$PROBES2Y{$probestring} = $ystring;
				$PROBES2SEQUENCE{$probestring} = $sequencestring;
			}
		}
	}
	print "finished querying the database\n";
}


sub writeTest{
	print "just writing test file $cdffile\n";
	open(OUT, "> $cdffile") or die "can't open output file $cdffile!";
	my @x;
	my @y;
	my @sequence;
	foreach my $key (keys %PROBES2TRANSCRIPT){
		@x = split /;/,$PROBES2X{ $key };
		print OUT $PROBES2TRANSCRIPT{ $key }."\t".$key."\t".(scalar(@x))."\n";
	}
	close(OUT);
}


######
## what I need:
## transcript_id,
## transcript_strand,
## exon_id,
## probe_chrom_start,
## probe_chrom_end
## get sum of nts of all exons downstream of exon
## example for forward: ENST00000375549, probe id 541251, exon ENSE00003609552, probe start: 111958630, probe end: 111958654 should be 1133+43=1176
## reverse ENST00000424047, probe id 499023, exon ENSE00001611854, 86785, 86809: singly exon.
##         ENST00000318333, 984421, ENSE00001258214, 93567099, 93567123: should be 803
sub getDistance3Prime{
  my @parms = @_;
  my $transcript_id=$parms[ 0 ];
  my $strand=$parms[ 1 ];
  my $exon_id=$parms[ 2 ];
  my $probe_chrom_start=$parms[ 3 ];
  my $probe_chrom_end=$parms[ 4 ];
  my $distance=0;
  $prep_query_get_exon->execute( $exon_id ) or die "can not execute query".$prep_query_get_exon->errstr();
  my (@exonstartend)=$prep_query_get_exon->fetchrow_array;
  if( $strand > 0 ){
    ## get the length of all exons downstream of the actual exon targeted by the probe.
    $prep_query_length_fw -> execute( $transcript_id, $probe_chrom_start ) or die "can not execute query ".$prep_query_length_fw -> errstr;
    my (@dummyresult) = $prep_query_length_fw->fetchrow_array;
    my $dist=0;
    ## if the probe is mapping the last exon it will return Null, obviously...
    if( @dummyresult ){
      if( defined( $dummyresult[ 0 ] ) ){
	$dist=$dummyresult[ 0 ];
      }
    }
    ## add the distance probe end to exon end to that.
    $distance= 0 + $dist + $exonstartend[ 1 ]-$probe_chrom_start;
  }else{
    ## get the length of all exons downstream of the actual exon targeted by the probe.
    $prep_query_length_rv -> execute( $transcript_id, $probe_chrom_start ) or die "can not execute query ".$prep_query_length_rv -> errstr;
    my (@dummyresult) = $prep_query_length_rv->fetchrow_array;
    my $dist=0;
    ## if the probe is mapping the last exon it will return Null, obviously...
    if( @dummyresult ){
      if( defined( $dummyresult[ 0 ] ) ){
	$dist=$dummyresult[ 0 ];
      }
    }
    ## add the distance probe end to exon end to that.
    $distance= 0 + $dist + $probe_chrom_end - $exonstartend[ 0 ];
  }
  return( $distance );
}


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


