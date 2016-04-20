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
use strict;
use Getopt::Std;

my %option=();
my $ensembl_database="core";
my $species="human";
my $transcript_id_column="transcript_id";
my $input_file;
my $output_file;
getopts("s:t:i:o:h",\%option);
if( defined( $option{ h } ) ){
  print( "\nannotateFileTranscript.pl version ".$script_version_.".\n" );
  print( "Retrieve various annotations for Ensembl transcript IDs provided in the input file using the Ensembl Perl API. Multiple values for a gene will be pasted into a single value, separated by a semicolon (;).\n\n");
  print( "Usage: perl annotateFileTranscript.pl -i:o:[t:hs:]\n\n" );
  print( "Parameters:\n" );
  print( "-i (required): the input file, i.e. a tabulator delimited text file that contains the Ensembl transcript IDs, one ID per row.\n" );
  print( "-o (required): the output file to which the annotation should be written. The output file will contain all columns from the input file with the annotations in additional columns.\n" );
  print( "-t (optional): the name of the column in the input file containing the Ensembl transcript IDs (defaults to transcript_id).\n" );
  print( "-h : print this help and exit\n" );
  print( "-s (optional): the species; defaults to human.\n" );
  exit 0;
}
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

my $gene_adaptor = $registry->get_adaptor( $species, $ensembl_database, "gene" );
#my $transcript_adaptor = $registry->get_adaptor( $species, $ensembl_database, "transcript" );

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
    @line = split /\t/,$_;
    if( $header==0 ){
      ## trying to find the column with the Ensembl gene ids in the input file.
      my $colcounter=0;
      foreach my $colname ( @line ){
	if( $colname eq $transcript_id_column ){
	  $transcript_id_idx = $colcounter;
	}
	$colcounter=$colcounter+1;
      }
      if( $transcript_id_idx < 0 ){
	die( "column $transcript_id_column not found in header of input file!\n" );
      }
      $header=1;
      print OUT $_."\tsymbol\tentrezgene\texternal_gene_id\tgene_biotype\tchromosome_name\tchromosome_strand\n";
    }
    else{
      my %symbol_hash=();
      my %entrezgene_hash=();
      my $symbol="";
      my $entrezgene="";
      my $external_name="";
      my $gene_biotype="";
      my $chromosome_name="";
      my $chromosome_strand="";
#      my $gene = $gene_adaptor->fetch_by_stable_id( $line[ $gene_id_idx ] );
      my $gene = $gene_adaptor->fetch_by_transcript_stable_id( $line[ $transcript_id_idx ] );
#      print_DBEntries( $gene->get_all_DBLinks( "EntrezGene" ) );
      if( defined( $gene ) ){
	my $all_entries = $gene->get_all_DBLinks( "EntrezGene" );
	foreach my $dbe ( @{$all_entries} ){
	  $symbol_hash{ $dbe->display_id } = "something";
	  $entrezgene_hash{ $dbe->primary_id } = "something";
	}
	$symbol = join( ";", sort keys %symbol_hash );
	$entrezgene = join( ";", sort keys %entrezgene_hash );
	$external_name = $gene->external_name;
	$gene->transform("chromosome");
	$chromosome_name = $gene->slice->seq_region_name;
	$gene_biotype=$gene->biotype;
	$chromosome_strand=$gene->strand;
      }
      else{
	print "transcript $line[ $transcript_id_idx ] not found in database\n";
      }
      print OUT $_."\t".$symbol."\t".$entrezgene."\t".$external_name."\t".$gene_biotype."\t".$chromosome_name."\t".$chromosome_strand."\n";
    }
  }
}



close( IN );
close( OUT );

