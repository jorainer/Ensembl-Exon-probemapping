#######
## 15.02.2011
my $script_version_="0.0.3";
##### description
## just a very cool perl script to annotate a tabulatir delimited text file with some additional columns from an. this is similar/identical to annotateFile.pl,
## just that here we use the transcript_id instead of the gene_id
## ensembl database.
# version 0.0.2: 02.02.2011: add chromosome name, strand and gene biotype
# version 0.0.3: 14.02.2011: use join instead of "concatenate"
# paramters:
# s: species, defaults to human.
# t: column in the input file containing the ensembl transcript IDs, defaults to transcript_id

#use lib '/home500/ensembl/temp/ensembl/ensembl/modules';
#use lib '/home500/ensembl/temp/58_37c/ensembl/ensembl/modules';
#use lib '/home/bioinfo/ensembl/67/API/ensembl/modules';
use strict;
use List::Util qw( min max );
use Getopt::Std;

my %option=();
my $ensembl_database="core";
my $species="human";
my $transcript_id_column="transcript_id";
my $input_file;
my $output_file;
getopts("s:g:i:o:",\%option);
if( defined( $option{ s } ) ){
  $species=$option{ s };
}
if( defined( $option{ t } ) ){
  $transcript_id_column=$option{ t };
}
if( defined( $option{ i } ) ){
  $input_file = $option{ i };
}
else{
  die( "no input file specified!\n" );
}

if( defined( $option{ o } ) ){
  $output_file = $option{ o };
}
else{
  die( "no output file specified!\n" );
}

use Bio::EnsEMBL::Registry;
my $user="anonuser";
my $host="madb.i-med.ac.at";
my $pass="";
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(-host => $host, -user => $user, -pass => $pass, -verbose => "1" );

#my $gene_adaptor = $registry->get_adaptor( $species, $ensembl_database, "gene" );
my $transcript_adaptor = $registry->get_adaptor( $species, $ensembl_database, "transcript" );

# Define a helper subroutine to print DBEntries
sub print_DBEntries
{
    my $db_entries = shift;

    foreach my $dbe ( @{$db_entries} ) {
        printf "\tXREF %s (%s) primary id %s\n", $dbe->display_id(), $dbe->dbname(), $dbe->primary_id;
    }
}

my $transcript_id_idx = -1;
my $header=0;
my @line;
open( IN, "< $input_file" ) or die "can't open file $input_file!\n";
open( OUT, "> $output_file" ) or die "can't open file $output_file for writing!\n";
while( <IN> ){
  chomp;
  if( /^#/ ){
  }
  else{
    ## get rid of quotes...
#    $_ =~ s/"//g;
    @line = split /\t/,$_;
    if( $header==0 ){
      ## trying to find the column with the Ensembl gene ids in the input file.
      my $colcounter=0;
      foreach my $colname ( @line ){
        $colname =~ s/"//g;
	if( $colname eq $transcript_id_column ){
	  $transcript_id_idx = $colcounter;
	}
	$colcounter=$colcounter+1;
      }
      if( $transcript_id_idx < 0 ){
	die( "column $transcript_id_column not found in header of input file!\n" );
      }
      $header=1;
      print OUT $_."\ttranscript_length\n";
    }
    else{
	## if we do have more than one transcript...
	my $length = 0;
	my @transcript_lengths=();
	my $transcript_line = $line[ $transcript_id_idx ];
	$transcript_line =~ s/"//g;
	my @transcripts = split /;/,$transcript_line;
	foreach my $transcript_id ( @transcripts ){
	    my $transcript = $transcript_adaptor->fetch_by_stable_id( $transcript_id );
	    if( defined( $transcript ) ){
		my $transcript_length = $transcript->length();
		push( @transcript_lengths, $transcript_length );
	    }
	}
	if( scalar (@transcript_lengths) > 0 ){
	    $length = min @transcript_lengths;
	}
#	else{
#	    print "transcript $line[ $transcript_id_idx ] not found in database\n";
#	}
	print OUT $_."\t".$length."\n";
    }
  }
}



close( IN );
close( OUT );

