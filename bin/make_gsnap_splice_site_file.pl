# This scripts generates a "known splice site file" that can be used for the GSNAP aligner.
# the process is the following:
# + get all genes from Ensembl for a species.
#   + foreach gene, get all transcripts.
#     + for each transcript, get all exons ($transcript->get_all_Exons, with the first exon being the most 5' exon), if more than one exon:
#       + if on + strand, use end, if on - strand use start.
my $script_version="0.0.1";

# loading additional modules
use DBI;
use strict;
use Getopt::Std;
use Bio::EnsEMBL::ApiVersion;
use Bio::EnsEMBL::Registry;

## default settings:
my $species = "human";
my $only_constitutive = 0;    ## evaluate only constitutive exons of a transcript.
## connect to the public Ensembl database.
my $user="anonymous";
my $host="ensembldb.ensembl.org";
my $mysql_port="5306";
my $pass="";
##
####

# Read option keys
my %option=();
getopts("hs:",\%option);
if( $option{ h } ){
    print( "\make_gsnap_splice_site_file.pl ".$script_version."\n" );
    print( "Create a known splice site file that can be used by the GSNAP aligner for all genes defined in Ensembl.\n\n" );
    print( "usage: make_gsnap_splice_site_file.pl -chs:\n" );
    print( "parameters:\n" );
    print( "-c if only constitutive exons of a transcript should be considered.\n" );
    print( "-h print this help.\n" );
    print( "-s (optional): the species. Defaults to human.\n\n" );
    exit;
}
if( defined( $option{ s } ) ){
  $species=$option{ s };
}
if( $option{ c } ){
  $only_constitutive = 1;
}

## connecting to the Ensembl database.
my $registry='Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(-host => $host, -user => $user, -pass => $pass, -port => $mysql_port, -verbose => "1" );
my $ensembl_version=$registry->software_version;
my $gene_adaptor = $registry->get_adaptor( $species ,"core","Gene");


## let's get rumbling
my $info_string="# make_gsnap_splice_site_file.pl version $script_version querying Ensembl version ".$ensembl_version."; ".localtime()."\n# Species: $species.\n# Only constitutive exons: $only_constitutive\n";
my $out_file = "SpliceSiteFile-Ensembl-".$ensembl_version."-".$species.".txt";
print "\n\n$info_string\nWriting output to:".$out_file."\n\n";
open( OUT, "> $out_file" ) or die "can't open file $out_file for writing!\n";
print OUT $info_string;
#do_test();
do_all();
close( OUT );
$registry -> disconnect_all();


## Takes a gene_id as input and processes that gene.
sub do_gene{
  my $current_gene_id = $_[ 0 ];
  my $orig_gene = $gene_adaptor->fetch_by_stable_id( $current_gene_id );
  my $do_transform = 1;
  my $gene = $orig_gene->transform( "chromosome" );
  if( !defined( $gene ) ){
    ## ach, gene is not on a chromosome...
    $gene = $orig_gene;
    $do_transform=0;
  }
  my $coord_system = $gene->coord_system_name;
  my $seq_name = $gene->seq_region_name;
  my $seq_strand = $gene->seq_region_strand;
  my @transcripts = @{ $gene->get_all_Transcripts };
  ## loop through the transcripts.
  foreach my $transcript ( @transcripts ){
    ## transform to chromsomal coordinates... just to be really sure...
    if( $do_transform==1 ){
      $transcript = $transcript->transform( "chromosome" );
    }
    my $transcript_id=$transcript->stable_id;
    ## get all exons.
    my @exons = @{ $transcript->get_all_Exons() };
    if( scalar( @exons ) > 1 ){
      ## obviously we can only define a splice site if the transcript has more than on exon...
      my $exon_idx = 1;
      my $donor_in_exon;      ## that's the nt position of the donor within the exon, will be exon_end for + and exon_start for - strand.
      my $acceptor_in_exon;   ## that's the nt position of the acceptor within the exon, will be exon_start for + and exon_end for - strand.
      my $donor_exon_id;
      my $acceptor_exon_id;
      my $intron_length;
      my $donor_string;
      my $acceptor_string;
      foreach my $exon (@exons){
	if( $do_transform==1 ){
	  $exon = $exon->transform( "chromosome" );
	}
	if( $exon_idx > 1 ){
	  ## ok, have defined a donor site, and will define the according acceptor site and write both to the output file.
	  if( $seq_strand < 0 ){
	    $acceptor_in_exon=$exon->end;
	    $intron_length=$donor_in_exon - $acceptor_in_exon - 1;
	    $acceptor_string="".($acceptor_in_exon+1)."..".$acceptor_in_exon;
	    $donor_string="".$donor_in_exon."..".($donor_in_exon-1);
	  }else{
	    $acceptor_in_exon=$exon->start;
	    $intron_length=$acceptor_in_exon - $donor_in_exon - 1;
	    $acceptor_string="".($acceptor_in_exon-1)."..".($acceptor_in_exon);  ## 1st nt of exon and last of intron
	    $donor_string="".$donor_in_exon."..".($donor_in_exon+1);           ## last of exon and first of intron
	  }
	  print OUT ">".$transcript_id.".".$current_gene_id.".".$donor_exon_id." ".$seq_name.":".$donor_string." donor ".$intron_length."\n";
	  print OUT ">".$transcript_id.".".$current_gene_id.".".$exon->stable_id." ".$seq_name.":".$acceptor_string." acceptor ".$intron_length."\n";
	}
	## define a donor site.
	if( $seq_strand < 0 ){
	  $donor_in_exon=$exon->start;
	}else{
	  $donor_in_exon=$exon->end;
	}
	$donor_exon_id = $exon->stable_id;
	$exon_idx++;
      }
    }
  }
}


####
## just to evaluate whether we get the same than in the README of gsnap.
## running the stuff for the genes ERBB2 and ERG.
## well, works, if coords from README are lifted to GRCH37.
## the NM_004448 transcript of ERBB2 is ENST00000540147.
## the NM_004449 transcript of ERG is ENST00000442448.
## while the example in the README seems to be wrong, the code above seems to work
## nicely (checked with the Ensembl genome browser).
sub do_test{
  my @genes = ( "ENSG00000141736", "ENSG00000157554" );
  foreach my $gene (@genes){
    do_gene( $gene );
  }
}


####
## defining splice junction for all genes defined in Ensembl.
sub do_all{
  my @genes = @{$gene_adaptor->list_stable_ids()};
  my $counta = 0;
  foreach my $gene (@genes){
    if( $counta % 1000 == 0 ){
      print localtime()."> Processed $counta of ".scalar(@genes)." genes.\n";
    }
    do_gene( $gene );
    $counta ++;
  }
}

