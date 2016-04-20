#!/usr/bin/perl
#####################################
my $script_version = "0.3.0";
##
# main script to start (cdna) alignments, i.e. to align probe sequences to cDNA.
#
# What are we doing here:
# 1) Create annotation database and initialize.
# 2) Align the probe sequences to the cDNA and ncrna fasta files.
# 3) Insert all alignments (up to xx missmatches; parameter number_mismatches_dbinsert) into the database.
# 4) Align all probe to the genome.
# 5) Re-evaluate all alignments and annotate them:
#    - for seq_type cdna: use the seq_id (transcript id) to get annotations from Ensembl.
#    - for genome... don't do anything, since we do not expect anything here...
# 6) Re-read the probe fasta file and summarize probe alignment informations.

#####################################
use IO::File;
use DBI;
use Getopt::Std;
use strict;
use warnings;
use Config::Simple;
use Parallel::ForkManager;
use Bio::EnsEMBL::ApiVersion;
use Bio::EnsEMBL::Registry;
use CustomCDF::DBAdaptor;
use CustomCDF::GeneAdaptor;
use CustomCDF::ProbeAdaptor;


my %option=();
getopts("k:e:h",\%option);
if( $option{ h } ){
    print( "\ncdna_alignment.pl version ".$script_version."\n" );
    print( "Aligns probe sequences of a microarray  against all cDNAs, genomic sequences and generates an annotation database for the specific microarray.\n\n" );
    print( "usage: cdna_alignment.pl -k[eh]\n" );
    print( "parameters:\n" );
    print( "-k A configuration file with all required settings (database connection settings, sequence files, parameters etc).\n" );
    print( "-e The Ensembl version. Optional.\n" );
    print( "-h print this help.\n\n" );
    exit;
}

if(!defined $option{ k }){
	die "the config file has to be specified (using the -k parameter)!\n";
}

## settings
my $cfg = new Config::Simple( $option{ k } );
my $host=$cfg->param( "database.mysql_host" );
my $username=$cfg->param( "database.mysql_user" );
my $password=$cfg->param( "database.mysql_password" );
my $dbname=$cfg->param( "database.dbprefix" );
my $ensembl_user=$cfg->param( "database.ensembl_user" );
my $ensembl_password=$cfg->param( "database.ensembl_password" );
if( $ensembl_password eq "empty" ){
    $ensembl_password="";
}
my $ensembl_host=$cfg->param( "database.ensembl_host" );
my $ensembl_database="core";

## general
my $exonerate=$cfg->param( "general.exonerate_binary" );
my $exonerate_out_path=$cfg->param( "general.exonerate_out_path" );
my $exonerate_out_file;
my $query_fasta_file=$cfg->param( "general.query_fasta_file" );
my $target_fasta_base_dir=$cfg->param( "general.target_fasta_base_dir" );
my $species=$cfg->param( "general.species" );
my $nr_row=$cfg->param( "general.nr_row" );

## alignment
my $number_mismatches_dbinsert=$cfg->param( "alignment.number_mismatches_dbinsert" );
my $min_align=$cfg->param( "alignment.min_alignment" );
my $probe_length=$cfg->param( "alignment.probe_length" );
my $number_mismatches_allowed=$cfg->param( "alignment.number_mismatches_allowed" );
my $chip_type=$cfg->param( "alignment.chip_type" );
my $max_number_GC=$cfg->param( "alignment.max_number_GC" );
my $max_number_genome_alignments=$cfg->param( "alignment.max_number_genome_alignments" );

## get the current Ensembl release from the file...
my $ensembl_version;
if( !defined $option{ e } ){
    my $cfg_file = $cfg->param( "general.config_file" );
    open( IN, "< $cfg_file" ) or die "can't open file $cfg_file!\n";
    while( <IN> ){
	chomp;
	if( /^ensembl/ ){
	    my @line = split /=/,$_;
	    $ensembl_version=$line[1];
	}
    }
    close( IN );
}else{
    $ensembl_version = $option{ e };
}
## check Ensembl version...
my $api_version="".software_version()."";
if( $ensembl_version ne $api_version ){
    die "The submitted Ensembl version (".$ensembl_version.") does not match the version of the Ensembl API (".$api_version."). Please configure your PERL5LIB environment variable to point to the correct API.";
}


## define the database name, output file name, path to the cDNA and genomic fasta files.
$dbname=$dbname."cdna_".$ensembl_version;
if ( -d $exonerate_out_path ){
}else{
    mkdir $exonerate_out_path;
}
$exonerate_out_file=$exonerate_out_path.$dbname.".out";
my $genomic_fasta_dir=$target_fasta_base_dir."/".$ensembl_version."/fasta/".$species."/dna/";
my $cdna_fasta_dir=$target_fasta_base_dir."/".$ensembl_version."/fasta/".$species."/cdna/";
my $ncrna_fasta_dir=$target_fasta_base_dir."/".$ensembl_version."/fasta/".$species."/ncrna/";

##################################### MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
print "This is cdna_alignment.pl version $script_version.\nUsing configuration file: $option{k}\n";

createDatabase();
createDBTables();
cdna_align_exonerate();
cdna_parse_exonerate();
ncrna_align_exonerate();
cdna_parse_exonerate();
align_all_genome();
create_index_probe_alignments(); ## finished with that table...
annotate_alignments();  ## looks good.
create_index_probe_to_exon_transcript();
create_exon_transcript_table();
probe_summary(); #OK
## define_nonoverlapping_regions...
define_exon_regions();
create_index_on_exon_region_tables();

print "\nFinished! Have fun!\n";
##
################

###################################### SUBS #####################################
## create database.
sub createDatabase{
    print localtime()."> Creating database $dbname...";
    my $dbcon = DBI->connect( "DBI:mysql:database=mysql;host=$host:user=$username:password=$password" ) or die "unable to connect to database!";
    my $dummyquery = $dbcon -> prepare( "CREATE database $dbname" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    $dbcon->disconnect();
    print "done.\n";
}

## create database tables.
sub createDBTables{
    my $dbcon = DBI->connect( "DBI:mysql:database=$dbname;host=$host:user=$username:password=$password" ) or die "unable to connect to database $dbname!";
    print localtime()."> Initializing database $dbname...";
    ## the core: the probe_alignments table... a little different here...
    my $dummyquery=$dbcon->prepare( "CREATE TABLE probe_alignments (probes_pk SERIAL,probe_id CHAR(7), x numeric, y numeric, seq_name CHAR(20), seq_type CHAR(10), start bigint, end bigint, seq_strand CHAR(1),missmatches numeric);" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    ## probe_to_exon_transcript v2
    $dummyquery=$dbcon->prepare( "CREATE TABLE probe_to_exon_transcript (transcript_id CHAR(20), gene_id CHAR(20), probes_fk bigint,chromosome_strand CHAR(1), probe_id CHAR(7), exon_id CHAR(25), probe_alignment_exon_start int, probe_alignment_exon_end int, position_of_alignment tinyint, is_junction_alignment tinyint);" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    print "done.\n";
    $dbcon->disconnect();
}

## start alignment of probes.
sub cdna_align_exonerate{
    print localtime()."> Starting alignment of probes in $query_fasta_file";
    my $dir=$cdna_fasta_dir;
    my $cdna_fasta_file="none";
    opendir (DIR, $dir) or die $!;
    while ( my $file=readdir( DIR ) ){
	if( $file=~m/cdna.all.fa$/ ){
	    $cdna_fasta_file=$file;
	}
    }
    closedir(DIR);
    my $targetfile = $cdna_fasta_dir.$cdna_fasta_file;
    print " against $cdna_fasta_file...";
    if(system("$exonerate -q $query_fasta_file -t $targetfile --querytype dna --targettype dna --forcescan t --showtargetgff no --showquerygff no --showvulgar no --showalignment no --showcigar no  --showsugar no --percent 88 --ryo \"sugar: \%S \%em\n\" > $exonerate_out_file")!=0){
        die "ERROR: exonerate call failed, last message: $! \n";
    }
    print "done.\n";
}


## parse results.
sub cdna_parse_exonerate{
    print localtime()."> Parsing $exonerate_out_file and storing its results to the database...";
    parse_exonerate( $exonerate_out_file, "cdna" );
    print "done.\n";
}

## start alignment of probes.
sub ncrna_align_exonerate{
    print localtime()."> Starting alignment of probes in $query_fasta_file";
    my $dir=$ncrna_fasta_dir;
    my $ncrna_fasta_file="none";
    opendir (DIR, $dir) or die $!;
    while ( my $file=readdir( DIR ) ){
	if( $file=~m/ncrna.fa$/ ){
	    $ncrna_fasta_file=$file;
	}
    }
    closedir(DIR);
    my $targetfile = $ncrna_fasta_dir.$ncrna_fasta_file;
    print " against $ncrna_fasta_file...";
    if(system("$exonerate -q $query_fasta_file -t $targetfile --querytype dna --targettype dna --forcescan t --showtargetgff no --showquerygff no --showvulgar no --showalignment no --showcigar no  --showsugar no --percent 88 --ryo \"sugar: \%S \%em\n\" > $exonerate_out_file")!=0){
        die "ERROR: exonerate call failed, last message: $! \n";
    }
    print "done.\n";
}


## parse results.
# sub cdna_parse_exonerate{
#     print localtime()."> Parsing $exonerate_out_file and storing its results to the database...";
#     parse_exonerate( $exonerate_out_file, "cdna" );
#     print "done.\n";
# }


## requires:
## exonerate output file
## seq_type
sub parse_exonerate{
    my @params = @_;
    my $exonerate_result = $params[ 0 ];
    my $seq_type = $params[ 1 ];
##########################################
## ok, what do i expect from the exonerate result?
## example line:
## sugar: probe:HuEx-1_0-st-v2:4009796;835:1566; 0 24 + X 112883320 112883344 + 111 1
## note: probe definition for 3' arrays are different: >probe:HG-U133_Plus_2:1007_s_at:718:317; , so there is no
##            ; to split the probe id etc!
## means: probe maps with nt 0 to 24 (24nt) to chromosome X, + strand from ... to ... with score 111 and 1 missmatch.
## with splitting at whitespaces i get:
## 0 sugar:
## 1 probe id
## 2 query start
## 3 query stop
## 4 query strand
## 5 target name
## 6 target start
## 7 target end
## 8 target direction
## 9 score
## 10 nr of missmatches
## hits to the + strand will result in : query start smaller then query end, query direction: + and the same for target
## to - strand it results in query start bigger then query end and query direction -, for the target strand:  start smaller end and + direction. (also + - is possible)
#my $count_length = 0;
    my $private_dbcon = DBI->connect( "DBI:mysql:database=$dbname;host=$host:user=$username:password=$password" ) or die "unable to connect to database $dbname!";
    my $prep_query=$private_dbcon -> prepare_cached("INSERT INTO probe_alignments (probe_id, x, y, seq_name, seq_type, start, end, seq_strand, missmatches) VALUES (?,?,?,?,?,?,?,?,?);");
    my $start;
    my $end;
    my $strand;
    my $xind;
    my $yind;
    my $probeid;
    my $missmatches;
    my $min_align_neg = -$min_align;
    ## read all lines from IN and write it to the database, if...
    open(IN, "$exonerate_result") or die "file $exonerate_result cannot be opened!\n";
    while(<IN>){
        chomp;
        if( /^sugar/){
            my @dummy = split /\s/;
	    my @tmp = split( /:/ , $dummy[ 1 ] );      # gives us probe HuEx-1_o... 400etc;835
	    my @tmp2 = split( /;/, $tmp[ 2 ] );
	    my @tmp3 = split( /-/, $tmp2[ 0 ] );
	    my $insertit=0;
	    ## treat ST and older 3' arrays differently; apparently Affymetrix changed the style of their FASTA files.
	    if( scalar( @tmp )==4 ){
		## got Exon or GeneST arrays
		$xind=$tmp2[ 1 ];
		$yind=$tmp[ 3 ];
		$yind =~ s/;//;
		$probeid=$tmp3[0];
	    }
	    if( scalar( @tmp )==5 ){
		## that's the old style.
		$xind=$tmp[ 3 ];
		$yind=$tmp[ 4 ];
		$yind =~ s/;//;
		## here I have to calculate the probe id...
		$probeid=1 + $nr_row * $yind + $xind;
	    }
	    my $diff = $dummy[ 3 ] - $dummy[ 2 ];
	    my $diff_neg = $dummy[ 2 ] - $dummy[ 3 ];
	    ## ok, now check if the hit fullfills the criterias to be included into the database.
	    if( $dummy[ 4 ] eq "+"  and $dummy[ 8 ] eq "+" and $dummy[ 10 ] <=$number_mismatches_dbinsert and $diff >=$min_align){
		#probe matches forward strand.
		$missmatches= $dummy[10] + $probe_length - $diff;	## in the case that exonerate has mapped less then the full probe sequence length (e.g. 24 instead of 25) we add the missing nucleotides to the missmatches.
		$start = $dummy[ 6 ]+1;		# have to add 1 to the start position, exonerate starts counting at 0 not at 1.
		$end = $dummy[ 7 ];
		$strand="+";
		if($missmatches <= $number_mismatches_dbinsert){
		    # just in case....
		    $insertit=1;
		}
	    }
	    if( $dummy[ 4 ] eq "-"  and $dummy[ 8 ] eq "+" and $dummy[ 10 ] <=$number_mismatches_dbinsert and $diff_neg >=$min_align){
		#probe matches reverse strand. query was reverse complemented
		$missmatches= $dummy[10] + $probe_length + $diff;	## in the case that exonerate has mapped less then the full probe sequence length (e.g. 24 instead of 25) we add the missing nucleotides to the missmatches.
		$start = $dummy[ 6 ]+1;
		$end = $dummy[ 7 ];
		$strand="-";
		if($missmatches <= $number_mismatches_dbinsert){
		    # just in case....
		    $insertit=1;
		}
	    }
	    if( $dummy[ 4 ] eq "+"  and $dummy[ 8 ] eq "-" and $dummy[ 10 ] <=$number_mismatches_dbinsert and $diff >= $min_align and $dummy[ 6 ] > $dummy[ 7 ] ){
		#probe matches reverse strand. target was reverse complemented by exonerate
		$missmatches= $dummy[10] + $probe_length - $diff;	## in the case that exonerate has mapped less then the full probe sequence length (e.g. 24 instead of 25) we add the missing nucleotides to the missmatches.
		$start = $dummy[ 7 ]+1;
		$end = $dummy[ 6 ];
		$strand="-";
		if($missmatches <= $number_mismatches_dbinsert){
		    # just in case....
		    $insertit=1;
		}
	    }
	    ## insert the exonerate result to the database.
	    if( $insertit==1 ){
		$prep_query -> execute( $probeid, $xind, $yind, $dummy[ 5 ], $seq_type, $start, $end, $strand, $missmatches ) or die "can not execute query ".$prep_query -> errstr;
	    }
        }
    }
    close(IN);
    $prep_query->finish();
    $private_dbcon->disconnect();
}


## align all probes also to the genome.
sub align_all_genome{
# Get a list of all genomic fasta files.
# start several processes to align against genome.
  my @genomic_fasta_files=();
  my $dir = $genomic_fasta_dir;
  opendir (DIR, $dir) or die $!;
  while ( my $file=readdir( DIR ) ){
    if( $file=~m/.fa$/ ){
      my @elements = split /\./,$file;
      ## fix for Ensembl >= 76: file name does no longer contain the Ensembl version, is thus shorter...
      if( $ensembl_version > 75 ){
	if( $elements[ 2 ] eq "dna" & $elements[ 3 ] eq "chromosome" & $elements[ 4 ]=~m/^(\d*|X|Y)$/ ){
	  push( @genomic_fasta_files, $file );
	}
      }else{
	if( $elements[ 3 ] eq "dna" & $elements[ 4 ] eq "chromosome" & $elements[ 5 ]=~m/^(\d*|X|Y)$/ ){
	  push( @genomic_fasta_files, $file );
	}
      }
    }
  }
  closedir(DIR);
  my $pm = Parallel::ForkManager->new( $cfg->param( "general.nrjobs" ) );
  ## ok, now align the probe sequences to these fasta files...
  foreach my $file (@genomic_fasta_files){
    my $pid = $pm->start and next;
    chromosome_align_exonerate( $file );  ## run the alignment of the probes.
    chromosome_parse_exonerate( $file );  ## insert the alignments to the database.
    $pm->finish;
  }
}

## align the probes to a chromosome
sub chromosome_align_exonerate{
    my @params=@_;
    my $outfile = $exonerate_out_path.$params[ 0 ];
    my $targetfile = $genomic_fasta_dir.$params[ 0 ];
    print localtime()."> Starting alignment of probes against $params[ 0 ].\n";
    if(system("$exonerate -q $query_fasta_file -t $targetfile --querytype dna --targettype dna --forcescan t --showtargetgff no --showquerygff no --showvulgar no --showalignment no --showcigar no  --showsugar no --percent 88 --ryo \"sugar: \%S \%em\n\" > $outfile")!=0){
        die "ERROR: exonerate call failed, last message: $! \n";
    }
    print localtime()."> Finished alignment of probes against $params[ 0 ].\n";
}

## parse the exonerate output and store alignments to the database.
sub chromosome_parse_exonerate{
    my @params=@_;
    my $exonerate_result = $exonerate_out_path.$params[ 0 ];
    print localtime()."> Parsing $exonerate_result and storing its results to the database.\n";
    parse_exonerate( $exonerate_result, "chromosome" );
    print localtime()."> Finished parsing $exonerate_result.\n";
}

## annotate alignments
sub annotate_alignments{
    ## we are not annotating all alignments, but only those with xx mismatches.
    ## also, we might consider to just annotate the cDNA alignments.
    ## we have to consider the alignment direction too!
    ## Workflow:
    ## 1) get all rows from probe_alignments against cdna on - strand (st array)
    ## 2) fetch transcript id and exons.
    ## 3) check where inside the exon the probe aligns.
    print localtime()."> Annotating all probe alignments against cDNA.\n";
    my $private_dbcon = DBI->connect( "DBI:mysql:database=$dbname;host=$host:user=$username:password=$password" ) or die "unable to connect to database $dbname!";
    ## connect to Ensembl
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_registry_from_db(
	-host => $ensembl_host,
	-user => $ensembl_user,
	-pass => $ensembl_password,
	-verbose => "1" );
    my $transcript_adaptor = $registry->get_adaptor( $species, 'Core','Transcript' );
    ## define the strand on which we expect the alignment.
    my $get_from_strand="+";
    if( $chip_type eq "sense" ){
	$get_from_strand="-";
    }
    my $insert_query_cached=$private_dbcon-> prepare( "insert into probe_to_exon_transcript (transcript_id, gene_id, probes_fk, chromosome_strand, probe_id, exon_id, probe_alignment_exon_start, probe_alignment_exon_end, position_of_alignment, is_junction_alignment) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);" ) or die $private_dbcon->errstr;
    ## fetch all alignments we want to annotate.
    my $fetch_alignments = $private_dbcon -> prepare( "select probes_pk, probe_id, seq_name, start, end from probe_alignments where seq_type='cdna' and seq_strand='".$get_from_strand."' and missmatches<=".$number_mismatches_allowed.";" ) or die $private_dbcon->errstr;
    $fetch_alignments->execute or die $private_dbcon->errstr();
    my $counter=0;
    while( my(@row) = $fetch_alignments->fetchrow_array ){
	## just some info...
	$counter++;
	if( ($counter % 5000) == 0 ){
	    print localtime()." > Processed $counter alignments\n";
	}
	my $probe_pk=$row[0];
	my $probe_id=$row[1];
	my $transcript_id=$row[2];
	my $alignment_start=$row[3];
	my $alignment_end=$row[4];
	my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);
	my $transcript_strand=$transcript->strand;
	my $gene = $transcript->get_Gene();
	my $gene_id = $gene->stable_id;
	my @exons = @{ $transcript->get_all_Exons( )};  #  -constitutive => 1 would return only constitutive exons...
	my $find_first=1;   ## this is 1 as long as we do not find the exon in which the probe alignment starts
	my $pos_of_alignment=0;
	for my $exon (@exons){
	    ## loop through exons and check if we do have the alignment within.
	    ## we do use only the exon ends, since the end of an exon is also the start of the next one.
	    my $exon_end=$exon->cdna_end( $transcript );
	    my $exon_start=$exon->cdna_start( $transcript );
	    my $probe_alignment_exon_start=0;
	    my $probe_alignment_exon_end=0;
	    if( $find_first==1 ){
		if( $alignment_start <= $exon_end ){
		    ## gotcha, alignment starts in this exon!
		    $find_first=0;
		    $pos_of_alignment=1;
		    $probe_alignment_exon_start=$alignment_start-$exon_start+1;
		    if( $alignment_end <= $exon_end ){
			## everything is within this exon.
			$probe_alignment_exon_end=$alignment_end-$exon_start+1;
			## > insert
			$insert_query_cached -> execute($transcript_id,
							$gene_id,
							$probe_pk,
							$transcript_strand,
							$probe_id,
							$exon->stable_id,
							$probe_alignment_exon_start,
							$probe_alignment_exon_end,
							$pos_of_alignment,0) or die $private_dbcon->errstr;
			last;  ## exit this loop.
		    }else{
			## got a junction-crossing probe, so the alignment end is in the next exon(s).
			$probe_alignment_exon_end=$exon_end-$exon_start+1;
			## > insert
			$insert_query_cached -> execute($transcript_id,
							$gene_id,
							$probe_pk,
							$transcript_strand,
							$probe_id,
							$exon->stable_id,
							$probe_alignment_exon_start,
							$probe_alignment_exon_end,
							$pos_of_alignment,1) or die $private_dbcon->errstr;
		    }
		}
	    }else{
		if( $alignment_end >= $exon_start ){ ## paranoia condition; make shure we have alignment covering that exon.
		    ## try to find the exon in which the alignment ends.
		    ## we will only get here if we do have a junction-crossing probe!
		    ## so, the start is in the previous exon, it is just the question whether the end is here
		    ## or if its also in the next.
		    $probe_alignment_exon_start=1;  ## alignment starts in previous exon, so in this exon we start at 1.
		    $pos_of_alignment=$pos_of_alignment+1;
		    if( $alignment_end <= $exon_end ){
			## that's it. insert and next.
			$probe_alignment_exon_end=$alignment_end-$exon_start+1;
			# > insert
			$insert_query_cached -> execute($transcript_id,
							$gene_id,
							$probe_pk,
							$transcript_strand,
							$probe_id,
							$exon->stable_id,
							$probe_alignment_exon_start,
							$probe_alignment_exon_end,
							$pos_of_alignment,1) or die $private_dbcon->errstr;
			last;  ## exit the loop.
		    }else{
			## still not the end... insert and move on.
			$probe_alignment_exon_end=$exon_end-$exon_start+1;  ## alignment ends at the exon end.
			# > insert
			$insert_query_cached -> execute($transcript_id,
							$gene_id,
							$probe_pk,
							$transcript_strand,
							$probe_id,
							$exon->stable_id,
							$probe_alignment_exon_start,
							$probe_alignment_exon_end,
							$pos_of_alignment,1) or die $private_dbcon->errstr;
		    }
		}
	    }
	}
    }
    ## finish all the stuff and cleanup
    $insert_query_cached->finish();
    $private_dbcon->disconnect();
    print localtime()."> Finished annotating all probe alignments against cDNA.\n";
}

## index on probe_alignments
## create index for columns:
## probe_id
## seq_name
## seq_type
## missmatches
sub create_index_probe_alignments{
    print localtime()."> Creating indices for table probe_alignments...";
    my $private_dbcon = DBI->connect( "DBI:mysql:database=$dbname;host=$host:user=$username:password=$password" ) or die "unable to connect to database $dbname!";
    my $dummyquery=$private_dbcon->prepare( "create index probe_alignments_probe_id_index ON probe_alignments (probe_id (7));" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    $dummyquery=$private_dbcon->prepare( "create index probe_alignments_seq_name_index ON probe_alignments (seq_name (16));" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    $dummyquery=$private_dbcon->prepare( "create index probe_alignments_seq_type_index ON probe_alignments (seq_type (4));" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    $dummyquery=$private_dbcon->prepare( "create index probe_alignments_missmatches_index ON probe_alignments (missmatches);" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    $private_dbcon->disconnect();
    print "done\n";
}

## index on probe_to_exon_transcript
## create index for columns:
## probe_id
## probes_fk
## transcript_id
## gene_id
sub create_index_probe_to_exon_transcript{
    print localtime()."> Creating indices for table probe_to_exon_transcript...";
    my $private_dbcon = DBI->connect( "DBI:mysql:database=$dbname;host=$host:user=$username:password=$password" ) or die "unable to connect to database $dbname!";
    my $dummyquery=$private_dbcon->prepare( "create index probe_to_exon_transcript_probe_id_index ON probe_to_exon_transcript (probe_id (7));" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    $dummyquery=$private_dbcon->prepare( "create index probe_to_exon_transcript_probes_fk_index ON probe_to_exon_transcript (probes_fk);" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    $dummyquery=$private_dbcon->prepare( "create index probe_to_exon_transcript_transcript_id_index ON probe_to_exon_transcript (transcript_id (15));" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    $dummyquery=$private_dbcon->prepare( "create index probe_to_exon_transcript_gene_id_index ON probe_to_exon_transcript (gene_id (15));" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    $dummyquery=$private_dbcon->prepare( "create index probe_to_exon_transcript_exon_id_index ON probe_to_exon_transcript (exon_id (15));" );
    $dummyquery->execute() or die "error: ".$dummyquery->errstr;
    $private_dbcon->disconnect();
    print "done\n";
}

#######
## this creates the "exon_transcript" table, each row representing one exon and its annotation to a transcript and gene
## the code was taken from the create_exon_transcript_table.pl
sub create_exon_transcript_table{
    print localtime()."> Create exon_transcript table.\n";
    my $private_dbcon = DBI->connect( "DBI:mysql:database=$dbname;host=$host:user=$username:password=$password" ) or die "unable to connect to database $dbname!";
    ## create the database table...
    my $private_query=$private_dbcon->prepare("create table exon_transcript ( exon_transcript_pk SERIAL, gene_id varchar(20), gene_name text, gene_biotype text, chromosome_name varchar(20), transcript_id varchar(20), transcript_chrom_start bigint, transcript_chrom_end bigint, transcript_chrom_strand tinyint, transcript_coding_chrom_start bigint, transcript_coding_chrom_end bigint, exon_id varchar(20), exon_chrom_start bigint, exon_chrom_end bigint, transcript_biotype text, coord_system text );") or die $private_dbcon->errstr();
    $private_query->execute() or die $private_dbcon->errstr;
    ## done
    ##
    my $insert_exon_transcript_cached = $private_dbcon->prepare_cached( "insert into exon_transcript (gene_id, gene_name, gene_biotype, chromosome_name, transcript_id, transcript_chrom_start, transcript_chrom_end, transcript_chrom_strand, transcript_coding_chrom_start, transcript_coding_chrom_end, exon_id, exon_chrom_start, exon_chrom_end, transcript_biotype, coord_system) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)" ) or die $private_dbcon->errstr();
    ##
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_registry_from_db(-host => $ensembl_host, -user => $ensembl_user, -pass => $ensembl_password, -verbose => "1" );
    my $gene_adaptor = $registry->get_adaptor( $species, $ensembl_database, "gene" );
    ## get all gene ids defined in the database...
    my @gene_ids = @{$gene_adaptor->list_stable_ids()};
    my $counta = 0;
    foreach my $gene_id ( @gene_ids ){
	$counta++;
	if( ($counta % 5000) == 0 ){
	    print localtime()."> Processed $counta genes\n";
	}
	my $orig_gene = $gene_adaptor->fetch_by_stable_id( $gene_id );
	if( defined $orig_gene ){
	    ## try to transform the gene to chromosome coord system, if not insert it anyway...
	    my $do_transform=1;
	    my $gene  = $orig_gene->transform("chromosome");
	    if( !defined $gene ){
	    	#next;
	    	## gene is not on known defined chromosomes!
	    	$gene = $orig_gene;
	    	$do_transform=0;
	    }
	    my $coord_system = $gene->coord_system_name;
	    my $chromosome_name = $gene->slice->seq_region_name;
	    my $gene_external_name= $gene->external_name;
	    my $gene_biotype = $gene->biotype;
	    my @transcripts = @{ $gene->get_all_Transcripts };
	    ## ok looping through the transcripts
	    foreach my $transcript ( @transcripts ){
		if( $do_transform==1 ){
		    $transcript = $transcript->transform("chromosome");	## just to be shure that we have the transcript in chromosomal coordinations.
		}
		my $transcript_start = $transcript->start;
		my $transcript_end = $transcript->end;
		my $transcript_strand = $transcript->strand;
		#my $transcript_chromosome= $transcript->seqname;
		my $transcript_cds_start = $transcript->coding_region_start;	## caution!!! will get undef if transcript is non-coding!
		if( !defined( $transcript_cds_start ) ){
		    $transcript_cds_start = "NULL";
		}
		my $transcript_cds_end = $transcript->coding_region_end;
		if( !defined( $transcript_cds_end ) ){
		    $transcript_cds_end = "NULL";
		}
		my $transcript_id = $transcript->stable_id;
		my $transcript_biotype = $transcript->biotype;
		my @exons = @{ $transcript->get_all_Exons() };
		foreach my $exon (@exons){
		    if( $do_transform==1 ){
			$exon->transform("chromosome");
		    }
		    my $exon_start = $exon->start;
		    my $exon_end = $exon->end;
		    my $exon_id = $exon->stable_id;
		    ## insert all this into the database...
		    ## use the prepared query.
		    $insert_exon_transcript_cached -> execute( $gene_id, $gene_external_name, $gene_biotype, $chromosome_name, $transcript_id, $transcript_start, $transcript_end, $transcript_strand, $transcript_cds_start, $transcript_cds_end, $exon_id, $exon_start, $exon_end, $transcript_biotype, $coord_system ) or die $private_dbcon->errstr;
		}
	    }
	}
    }
    $insert_exon_transcript_cached->finish;
    ## create indices...
    $private_query = $private_dbcon->prepare("create index exon_transcript_exon_id_idx on exon_transcript ( exon_id( 15 ) );") or die $private_dbcon->errstr();
    $private_query -> execute() or die $private_dbcon->errstr;
    $private_query = $private_dbcon->prepare("create index exon_transcript_gene_id_idx on exon_transcript ( gene_id( 15 ) );") or die $private_dbcon->errstr();
    $private_query -> execute() or die $private_dbcon->errstr;
    $private_query = $private_dbcon->prepare("create index exon_transcript_transcript_id_idx on exon_transcript ( transcript_id( 15 ) );") or die $private_dbcon->errstr();
    $private_query -> execute() or die $private_dbcon->errstr;
    ## done
    $private_dbcon->disconnect();
    print localtime()."> Done creating exon_transcript_table\n.";
}



#######
## this was the previous processProbemapping.pl
## we are summarizing all information for each probe.
## thus we re-read the fasta file and check for each probe how often it matches the genome etc.
## what do i expect form the fasta file?
## example:
## >probe:HuEx-1_0-st-v2:494998;917:193; ProbeSetID=....
## CAGGCCTTTTAATTTTTT
## so i will get the platform, probe_id, x and y coordinates from the first line, the sequence from the second
sub probe_summary{
    my $nr_chrom_map;
    my $nr_chrom_map_mm1;
    my $nr_chrom_map_mmall;
    my $gc_count;
    my $nr_gene_map;
    my $nr_gene_map_mm;  ## that's new: count the number of times a probe aligns to genes with one mismatches.
    my $probe_id;
    my $platform;
    my $xind;
    my $yind;
    my $sequence;
    print localtime()."> Gathering all alignment informations for each probe.\n";
    open(IN, "$query_fasta_file") or die "file $query_fasta_file cannot be opened!\n";
    my $private_dbcon = DBI->connect( "DBI:mysql:database=$dbname;host=$host:user=$username:password=$password" ) or die "unable to connect to database $dbname!";
    ## create table
    my $create_query = $private_dbcon->prepare("CREATE TABLE probe_information (probe_id CHAR(7), x smallint, y smallint, platform text, sequence text, nr_chrom_map mediumint, nr_chrom_map_mm1 mediumint, nr_chrom_map_mmall mediumint, nr_gene_map mediumint, nr_gene_map_mm mediumint, gc_count tinyint, cel_index mediumint);") or die $private_dbcon->errstr();
    $create_query->execute or die $private_dbcon->errstr;
    my $prep_insert_query = $private_dbcon->prepare_cached("INSERT INTO probe_information (probe_id, x, y, platform, sequence, nr_chrom_map, nr_chrom_map_mm1, nr_chrom_map_mmall, nr_gene_map, nr_gene_map_mm, gc_count, cel_index) VALUES (?,?,?,?,?,?,?,?,?,?,?,?);");
    my $pq_probe_matches_genome = $private_dbcon->prepare_cached( "select missmatches,count(*) from probe_alignments where probe_id=? and seq_type='chromosome' group by missmatches;" ) or die $private_dbcon->errstr();
    my $pq_probe_matches_gene = $private_dbcon->prepare_cached( "select count(distinct gene_id) from probe_to_exon_transcript where probe_id=?" ) or die $private_dbcon->errstr();
    my $pq_probe_matches_gene_with_mm = $private_dbcon->prepare_cached( "select count(distinct gene_id) from (select distinct seq_name from probe_alignments where probe_id=? and missmatches>0 and seq_type='cdna') as firstquery join exon_transcript on (transcript_id=seq_name)" ) or die $private_dbcon->errstr();
    ##
    my $counta=0;
    while(<IN>){
        chomp;
	$counta++;
	if( ($counta % 10000) == 0 ){
	    print localtime()."> Processed $counta probes\n";
	}
        if( /^>probe/){
	    ## defining the probe_id and all other required informations
	    my @splittedline = split /\s/;
	    my @toprocess = split( /:/, $splittedline[ 0 ] );
	    $platform = $toprocess[ 1 ];

	    my @processed = split( /;/, $toprocess[ 2 ] );
	    my @andagain = split( /-/, $processed[ 0 ] );
	    ## treat ST and older 3' arrays differently; apparently Affymetrix changed the style of their FASTA files.
	    if( scalar( @toprocess )==4 ){
		## got Exon or GeneST arrays
		$xind=$processed[ 1 ];
		$yind=$toprocess[ 3 ];
		$yind =~ s/;//;
		$probe_id=$andagain[0];
	    }
	    if( scalar( @toprocess )==5 ){
		## that's the old style.
		$xind=$toprocess[ 3 ];
		$yind=$toprocess[ 4 ];
		$yind =~ s/;//;
		## here I have to calculate the probe id...
		$probe_id=1 + $nr_row * $yind + $xind;
	    }
	}
	if( /^(A|T|C|G)/ ){
		## assuming we have now the sequence of the probe defined in the previous loop.
		$sequence=$_;
		## counting occurrences of G and C:
		$gc_count = 0;
		$gc_count = $gc_count + ($sequence =~ tr/G//);
		$gc_count = $gc_count + ($sequence =~ tr/C//);

		## initializing the variables.
		$nr_chrom_map=0;
		$nr_chrom_map_mm1=0;
		$nr_chrom_map_mmall=0;
		$nr_gene_map=0;
		$nr_gene_map_mm=0;
		## checking how often the probe maps to the genome.
		$pq_probe_matches_genome->execute($probe_id) or die $private_dbcon->errstr;
		while(my(@row) = $pq_probe_matches_genome->fetchrow_array) {
			## if missmatches=0 count to $nr_chrom_map, else to the mismatch counts...
			if($row[ 0 ] == 0){
				$nr_chrom_map = $nr_chrom_map + $row[ 1 ];
			}
			if( $row[ 0 ] == 1){
				## got counts for 1 missmatch.
				$nr_chrom_map_mm1=$nr_chrom_map_mm1 + $row[ 1 ];
				$nr_chrom_map_mmall = $nr_chrom_map_mmall + $row[ 1 ];
			}
			if( $row[ 0 ] > 1){
				## got more than one mismatch
				$nr_chrom_map_mmall = $nr_chrom_map_mmall + $row[ 1 ];
			}
		}
		## perfect alignment of probe to gene: use probe_to_exon_transcript table
		$pq_probe_matches_gene->execute($probe_id) or die $private_dbcon->errstr();
		while(my(@row) = $pq_probe_matches_gene->fetchrow_array) {
		    $nr_gene_map = $row[ 0 ];
		}
		## alignment with mismatches to gene: use join of probe_alignments and exon_transcript table
		#$private_query = $private_dbcon->prepare( "select count(distinct gene_id) from (select distinct seq_name from probe_alignments where probe_id=\'$probe_id\' and missmatches>0 and seq_type='cdna') as firstquery join exon_transcript on (transcript_id=seq_name)" ) or die $private_dbcon->errstr();
		$pq_probe_matches_gene_with_mm->execute($probe_id) or die $private_dbcon->errstr();
		while(my(@row) = $pq_probe_matches_gene_with_mm->fetchrow_array) {
		    $nr_gene_map_mm = $row[ 0 ];
		}

		## ok, now that i've gathered all required informations let's them insert into the database...
		$prep_insert_query -> execute( $probe_id,  $xind, $yind, $platform, $sequence, $nr_chrom_map, $nr_chrom_map_mm1, $nr_chrom_map_mmall, $nr_gene_map, $nr_gene_map_mm, $gc_count, getIndex( $xind, $yind ) ) or die "can not execute query ".$prep_insert_query -> errstr;
	}
    }
    close( IN );
    $prep_insert_query->finish();
    $pq_probe_matches_genome->finish;
    $pq_probe_matches_gene->finish;
    $pq_probe_matches_gene_with_mm->finish;
    ## create index.
    my $atlast = $private_dbcon->prepare( "create index probe_information_probe_id_idx on probe_information (probe_id(7))" ) or die $private_dbcon->errstr();
    $atlast->execute() or die $private_dbcon->errstr;
    $private_dbcon->disconnect();
    print localtime()."> Finished\n.";
}

sub getIndex{
	my @xy = @_;
	return ( ( $xy[ 0 ] + $xy[ 1 ] * $nr_row ) +1 );
}


############### SUBS for exon_region
## what are we going to do here?
# 1) create the database tables.
# 2) on a gene by gene basis do:
#    - get all exons of the gene and define the non-overlapping regions.
#    - get all probes matching the gene.
#    - check whether probe matches within exon-region: to do that I have to get
sub define_exon_regions{
  ## create the tables
  create_exon_region_tables();
  print localtime()."> Started creating non-overlapping exon regions for each gene and annotating probes to them.\n";
  my $db_adaptor = CustomCDF::DBAdaptor->new( host=>$host, username=>$username, password=>$password, dbname=>$dbname );
  ## prepare insert queries.
  my $prep_query_insert_regions_to_transcripts = $db_adaptor->dbcon->prepare_cached( "insert into ensembl_exon_regions_to_transcripts ( transcript_id, region_id, exon_id, inside_coding ) values ( ?, ?, ?, ? );" );
  ## insert exon_region
  my $prep_query_insert_region = $db_adaptor->dbcon->prepare_cached( "insert into ensembl_exon_regions ( gene_id, chromosome_name, chromosome_strand, region_id, region_chrom_start, region_chrom_end, exons, ensembl_evidence, astd_evidence, potentially_coding) values ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ? );" );
  ## insert probe exon region mapping
  my $prep_query_insert_probe_to_exon_regions = $db_adaptor->dbcon->prepare_cached( "insert into ensembl_probe_to_exon_regions ( probe_id, region_id, probe_chrom_start ) values ( ?, ?, ? );" );
  ### define query to fetch transcripts and coding start for an exon.
  my $prep_query_transcript_for_exon = $db_adaptor->dbcon->prepare_cached( "select transcript_id, transcript_coding_chrom_start, transcript_coding_chrom_end from exon_transcript where exon_id=?" );

  ## define the exon_regions for all genes.
  ## the point is whether it is at all possible to do...
  ## in "probe_to_exon_transcript" table I've got the "good" probes with their alignment
  ## relative within the exon. this in combination with exon start/end from "exon_transcript"
  ## could be used to assign probes to exon_regions.
  my $gene_query = $db_adaptor->dbcon->prepare( "select distinct gene_id from exon_transcript;" ) or die $db_adaptor->dbcon->errstr();
  $gene_query->execute() or die $db_adaptor->dbcon->errstr();
  my $gene_counter=0;
  my $gene_adaptor = $db_adaptor->get_gene_adaptor();
  my $probe_adaptor = $db_adaptor->get_probe_adaptor();
  while( my @row = $gene_query->fetchrow_array){
    ## things we are going to do here:
    my $gene_id= $row[ 0 ];
    my @STARTENDS = ();   ## array that will contain start and end coodinates of exons.
    my @WHAT = ();        ## array with 1 and -1 indicating whether the element at the same position in @STARTENDS is a start or end respectively.
    my @REGIONSTART= ();
    my @REGIONEND= ();
    my %EXON2START=();
    my %EXON2END=();
    my %PROBE2START=();
    my $chromosome_name;
    my $chromosome_strand;
    ## 1) get all exons for that gene (table "exon_transcript") and define non-overlapping
    ##    exon regions.
    my @exons = $gene_adaptor->fetch_exons_for_gene( gene_id=>$gene_id );  ## get exons ordered by chromosome start.
    $chromosome_name=$exons[ 0 ]->chromosome_name;
    $chromosome_strand=$exons[ 0 ]->strand;
    foreach my $exon (@exons){
      $EXON2START{ $exon->exon_id } = $exon->start;
      $EXON2END{ $exon->exon_id } = $exon->end;
      push( @STARTENDS, $exon->start );
      push( @STARTENDS, $exon->end );
      push( @WHAT, 1 );
      push( @WHAT, -1 );
    }

    #######
    # ok, have to get a ordering index for the start/end coordinates
    # which we will use then for the definition of the non-overlapping
    # exon regions.
    my @list_order = sort{ $STARTENDS[ $a ] <=> $STARTENDS[ $b ] } 0 .. $#STARTENDS;
    @STARTENDS = @STARTENDS[ @list_order ];
    @WHAT = @WHAT[ @list_order ];

    # that's fine, now i can loop through this array and define the
    # regions!
    my $the_counter=0;
    for my $i ( 0 .. $#STARTENDS ){
      if( $the_counter==0 ){
	push( @REGIONSTART, $STARTENDS[ $i ] );
	$the_counter=1;
      }
      else{
	push( @REGIONEND, $STARTENDS[ $i ] );
	# add +1 if we have a exon start, -1 if it is the end of an exon.
	$the_counter = $the_counter + $WHAT[ $i ];
	if( $the_counter > 0 ){
	  push( @REGIONSTART, $STARTENDS[ $i ] );
	}
    }
    }
    # non-overlapping regions are now defined!
    #######

    ## get all "good" probes for that gene.
    ## * targets at least one exon of the gene.
    ## * 1 or less genomic alignments
    ## * no genomic alignment with 1 mismatch.
    ## * less than xx G-Cs.
    ## Note: not using fetch_probes_for_gene(gene_id=>$gene_id, load_alignments=>1, seq_type=>"='chromosome'") because we can have some gene that are not on "conventional" chromosomes, but patches etc.
    my @probes = $probe_adaptor->fetch_probes_for_gene( gene_id=>$gene_id, load_genomic_alignments=>1 );
    ## now we have an array with Probe objects, each having one Alignment object (not more due to the restriction to a single genomic alignment).
    foreach my $probe (@probes){
      if( $probe->gc_count <= $max_number_GC ){
	my @alignments = $probe->alignments();
	my @starts = $alignments[ 0 ]->start();
	$PROBE2START{ $probe->id() } = $starts[ 0 ];
      }
    }

    ## todo:
    ## * have to define whether exon region is potentially coding.
    ## * do we have probes within the region?
    my $region_counter = 0;
    for my $i ( 0 .. $#REGIONSTART ){
      my $current_start = $REGIONSTART[ $i ];
      my $current_end = $REGIONEND[ $i ];
      my $has_probes = 0;
      my $potentially_coding = 0;
      my $exons = "";
      my $sep = "";
      my %PROBES = ();
      $region_counter++;
      my $region_id = $gene_id.":".$region_counter;
      my $ensembl_evidence=1;
      my $astd_evidence=0;
      foreach my $key ( keys %EXON2START ) {
	if( $EXON2START{ $key } <= $current_start && $EXON2END{ $key } >= $current_end ){
	  ## ok, exon contains region (region is within exon boundaries)
	  ## get all transcripts of this exon, check if the transcript is "coding",
	  ## if yes, check if the region is within the coding region of this transcript
	  ## and insert the mapping of the exon_region to the transcripts to the database.
	  $prep_query_transcript_for_exon -> execute( $key ) or die $prep_query_transcript_for_exon->errstr();
	  while( my @row = $prep_query_transcript_for_exon->fetchrow_array ){
	    ## row has: 0 transcript_id, 1 transcript_coding_chrom_start, 2 transcript_coding_chrom_end
	    my $inside_coding="no";
	    ## is the transcript coding??
	    if( defined( $row[ 1 ] ) ){
	      ## check if the region is fully within the coding regions
	      if( $current_start >= $row[ 1 ] && $current_end <= $row[ 2 ] ){
		$inside_coding="full";
		$potentially_coding=1;
	      }
	      elsif( ( $current_start <= $row[ 1 ] && $current_end >= $row[ 1 ] ) | ( $current_start <= $row[ 2 ] && $current_end >= $row[ 2 ] ) ){
		## region is partially within the coding region of the transcript
		$inside_coding="partially";
		$potentially_coding=1;
	      }
	    }
	    ## what remains is to insert this to the database.
	    $prep_query_insert_regions_to_transcripts->execute( $row[ 0 ], $region_id, $key, $inside_coding ) or die $prep_query_insert_regions_to_transcripts->errstr;
	  }
	  $exons = $exons.$sep.$key;
	  $sep=";";
	}
      }
      ## insert this friendly non-overlapping exon-region to the database...
      $prep_query_insert_region->execute( $gene_id, $chromosome_name, $chromosome_strand, $region_id, $current_start, $current_end, $exons, $ensembl_evidence, $astd_evidence, $potentially_coding ) or die $prep_query_insert_region->errstr;

      # search for probes within the region...
      foreach my $key ( keys %PROBE2START ){
	if( $PROBE2START{ $key } >= $current_start && ( $PROBE2START{ $key } + 24 ) <= $current_end ){
	  $has_probes = 1;
	  $PROBES{ $key } = $PROBE2START{ $key };
	  delete( $PROBE2START{ $key } );    # this will speed up the
	  # next loops...
	}
      }
      # ok, now insert the region with assigned probes... to the
      # database, but only if the region is targeted by probes
      if( $has_probes == 1 ){
	# what remains to do is to insert the probe to region mapping to the database.
	foreach my $key ( keys %PROBES ) {
	  $prep_query_insert_probe_to_exon_regions -> execute( $key, $region_id, $PROBES{ $key } ) or die $prep_query_insert_probe_to_exon_regions->errstr;
	}
      }
    }
    ## just some feedback...
    $gene_counter++;
    if( ( $gene_counter % 5000 ) == 0 ){
      print localtime()."> Processed ".$gene_counter." genes\n";
    }
  }
  print localtime()."> Done with the definition of non-overlapping exon regions.\n";
  $prep_query_insert_regions_to_transcripts->finish;
  $prep_query_insert_region->finish;
  $prep_query_insert_probe_to_exon_regions->finish;
  $prep_query_transcript_for_exon->finish;

#  $prep_query_exon_transcript->finish;
#  $private_dbcon->close;
}

## what tables are we creating?
## first: exon_regions: defines all non-overlapping exon_regions.
## second: exon_regions_to_transcripts: mapping between transcripts and exon regions.
## third: probe_to_exon_regions: probes assigned to exon regions.
sub create_exon_region_tables{
  print localtime()."> Creating exon region database tables...";
  my $private_dbcon = DBI->connect( "DBI:mysql:database=$dbname;host=$host:user=$username:password=$password" ) or die "unable to connect to database $dbname!";
  # create the exon_regions table
  my $query = $private_dbcon->prepare( "create table ensembl_exon_regions ( exon_regions_pk SERIAL, gene_id varchar(20), chromosome_strand tinyint, chromosome_name varchar(10), region_id varchar( 23 ), region_chrom_start bigint, region_chrom_end bigint, exons text, ensembl_evidence tinyint, astd_evidence tinyint, potentially_coding tinyint )" ) or die $private_dbcon->errstr();
  $query -> execute() or die $private_dbcon->errstr;
  # create the regions_to_transcripts table
  $query = $private_dbcon->prepare( "create table ensembl_exon_regions_to_transcripts ( region_id varchar(23), transcript_id varchar(20), exon_id varchar(20), inside_coding varchar(9) )" ) or die $private_dbcon->errstr();
  $query -> execute() or die $private_dbcon->errstr;
  # create the probe to exon regions table
  $query = $private_dbcon->prepare( "create table ensembl_probe_to_exon_regions ( probe_to_exon_regions_pk SERIAL, probe_id varchar(20), region_id varchar(23), probe_chrom_start bigint );" ) or die $private_dbcon->errstr();
  $query -> execute() or die $private_dbcon->errstr;
  $private_dbcon->disconnect();
  print "done\n";
}

sub create_index_on_exon_region_tables{
  print localtime()."> Creating index for exon region database tables...";
  my $private_dbcon = DBI->connect( "DBI:mysql:database=$dbname;host=$host:user=$username:password=$password" ) or die "unable to connect to database $dbname!";
  ## ensembl_exon_regions
  my $index_query = $private_dbcon->prepare( "create index exon_regions_gene_id_idx on ensembl_exon_regions( gene_id(20) );" ) or die $private_dbcon->errstr();
  $index_query->execute() or die $private_dbcon->errstr();
  $index_query = $private_dbcon->prepare( "create index exon_regions_region_id_idx on ensembl_exon_regions( region_id(20) );" ) or die $private_dbcon->errstr();
  $index_query->execute() or die $private_dbcon->errstr();

  ## ensembl_exon_regions_to_transcripts
  $index_query = $private_dbcon->prepare( "create index regions_to_transcripts_transcript_id on ensembl_exon_regions_to_transcripts(transcript_id(20));" ) or die $private_dbcon->errstr;
  $index_query -> execute() or die $private_dbcon->errstr;
  $index_query = $private_dbcon->prepare( "create index regions_to_transcripts_region_id on ensembl_exon_regions_to_transcripts(region_id(20));" ) or die $private_dbcon->errstr;
  $index_query -> execute() or die $private_dbcon->errstr;
    $index_query = $private_dbcon->prepare( "create index regions_to_transcripts_exon_id on ensembl_exon_regions_to_transcripts(exon_id(20));" ) or die $private_dbcon->errstr;
  $index_query -> execute() or die $private_dbcon->errstr;

  ## probe_to_exon_regions table: probe_id, region_id
  $index_query = $private_dbcon->prepare( "create index probe_to_exon_regions_probe_id on ensembl_probe_to_exon_regions(probe_id(20));" ) or die $private_dbcon->errstr;
  $index_query -> execute() or die $private_dbcon->errstr;
  $index_query = $private_dbcon->prepare( "create index probe_to_exon_regions_region_id on ensembl_probe_to_exon_regions(region_id(20));" ) or die $private_dbcon->errstr;
  $index_query -> execute() or die $private_dbcon->errstr;
  $private_dbcon->disconnect();
  print "done\n";
}





## try code with perl shell:
## psh

## annotate alignments to transcripts and genes.
## Notes:
# use Bio::EnsEMBL::Registry;
#
# Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
#  $transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor( 'Human', 'Core','Transcript' );
# $transcript = $transcript_adaptor->fetch_by_stable_id('ENST00000201961');
# $exon_adaptor = Bio::EnsEMBL::Registry->get_adaptor( 'Human', 'Core','Exon' );
# my @exons = @{ $transcript->get_all_Exons( )};  #  -constitutive => 1 would return only constitutive exons...
# for my $exon (@exons){
#    print( $exon->stable_id.": ".$exon->cdna_start($transcript)."-".$exon->cdna_end($transcript)."\n");
# }
##
# Transcript:
# might work: transcript get exons, on the exons call start, if not transformed to chromosome coordinate we might get
# relative coordinates to the transcript.

