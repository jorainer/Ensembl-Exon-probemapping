#!/usr/bin/perl
# script to fetch conservation constrained elements from Ensembl compara.
# options:
# -s: species.
# -m: method link type: e.g. EPO_LOW_COVERAGE, or EPO, GERP_CONSTRAINED_ELEMENT etc
# -c: species set name; collection of species used in the conservation analysis: e.g. amniotes, mammals, fish, primates.
use IO::File;
use DBI;
use Getopt::Std;
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

#use lib "/home/bioinfo/ensembl/72/API/ensembl/modules";            #  <- change here Ensembl API
#use lib "/home/bioinfo/ensembl/72/API/ensembl-compara/modules";            #  <- change here Ensembl API
my $script_version_="0.0.2";

## check Ensembl version...
use Bio::EnsEMBL::ApiVersion;
my $api_version="".software_version()."";

## Ensembl connection settings:
my $mysql_host="ensembldb.ensembl.org";
my $mysql_user="anonymous";
my $mysql_passwd="";
my $mysql_port="5306";
#my $mysql_host="madb.i-med.ac.at";
#my $mysql_user="anonuser";
#my $mysql_passwd="";
#my $mysql_port="3306";
##
### connect to the database
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(-host => $mysql_host, -user => $mysql_user, -pass => $mysql_passwd, -port => $mysql_port ,-verbose => "1" );


my $species="human";
my $species_set="";
my $method_link_type="";
my %option=();
getopts("s:m:c:",\%option);
if( defined( $option{ s } ) ){
  $species = $option{ s };
}
if( defined( $option{ m } ) ){
  $method_link_type = $option{ m };
}
else{
  print "No method link defined (option -m)! Options are:\n";
  printMethodLinks();
  exit;
}
if( defined( $option{ c } ) ){
  $species_set = $option{ c };
}
else{
  print "No species set defined (option -c)! Options are:\n";
  printSpeciesSetsByMethodLink();
  exit;
#  die "The species set has to be defined using the -c parameter. Options are: amniotes, birds, fish, primates, mammals.\n";
}

## compile the output file name:
my $out_file=$api_version."_".$species."_".$method_link_type."_".$species_set.".txt";

my $method_link_species_set_adaptor=$registry->get_adaptor( 'Multi', 'compara', 'MethodLinkSpeciesSet' );
my $method_link_species_set = $method_link_species_set_adaptor->fetch_by_method_link_type_species_set_name( $method_link_type, $species_set );
throw("Unable to find method_link_species_set") if (!defined($method_link_species_set));
my $cons_element_adaptor=$registry->get_adaptor( 'Multi', 'compara', 'ConstrainedElement' );

my $info_string="# This is fetch_constrained_elements_by_species.pl version $script_version_; settings:\n# species:\t$species\n# Ensembl version:\t$api_version\n# species_set:\t$species_set\n# method_link_type:\t$method_link_type\n# name:\t".$method_link_species_set->name()."\n";
print $info_string;
open( OUT, "> $out_file" ) or die "cannot open output file $out_file!\n";
print OUT $info_string;



my $slice_adaptor=$registry->get_adaptor( $species, "core", "slice" );

##############
## debugging...
#my $test_slice = $slice_adaptor->fetch_by_region("chromosome", "15", 76628758, 76635191 );
#my $test_ces = $cons_element_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice( $method_link_species_set, $test_slice );
#print "got".@$test_ces."\n";
##
##############

my $slices = $slice_adaptor->fetch_all( "chromosome" );  # get all chromosomes.
print OUT "chromosome\tstrand\tstart\tend\tscore\tp_value\n";
#foreach my $slice ( @{ $slice_adaptor->fetch_all( "chromosome" ) } ){
foreach my $slice (@$slices){
#  print "".$slice->name()."\n";
  my $slice_name=$slice->seq_region_name();
  print "processing chromosome ".$slice_name.": ";
  my $ce_listref=$cons_element_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice( $method_link_species_set, $slice );
  print "found ".@$ce_listref." elements\n";
  foreach my $constrained_element (@$ce_listref){
    my $strand = $constrained_element->strand();
    my $start = $constrained_element->start;
    my $end = $constrained_element->end();
    my $score = $constrained_element->score();
    my $pvalue = $constrained_element->p_value();
    if( !defined( $pvalue ) ){
      $pvalue="";
    }
    #my $tax = $constrained_element->taxonomic_level;
    print OUT $slice_name."\t".$strand."\t".$start."\t".$end."\t".$score."\t".$pvalue."\n";
  }
}


close( OUT );
$registry -> disconnect_all();

###############
## list all possible method links.
sub printMethodLinks{
  my $method_link_species_set_adaptor=$registry->get_adaptor( 'Multi', 'compara', 'MethodLinkSpeciesSet' );
  #my $species_set_adaptor=$registry->get_adaptor( 'Multi', 'compara', 'SpeciesSet' );

  my $ref_all = $method_link_species_set_adaptor->fetch_all();
#my $ref_all = $method_link_species_set_adaptor->fetch_all_by_method_link_type( $method_link_type );
  foreach my $mlssa (@$ref_all){
    print "name: ".$mlssa->name." method link type: ".$mlssa->method->type()."\n";
#    print "name: ".$mlssa->name." method link type: ".$mlssa->method_link_type."\n";  ## pre 70
#    print "" species set id:".$mlssa->species_set_id." source:".$mlssa->source."\n";
#    my $specset=$species_set_adaptor->fetch_by_dbID( $mlssa->species_set_id );
#    print "species set:".$specset->name()."\n\n";
  }
}


###############
## list all species sets for a method link
sub printSpeciesSetsByMethodLink{
  my $method_link_species_set_adaptor=$registry->get_adaptor( 'Multi', 'compara', 'MethodLinkSpeciesSet' );
  my $species_set_adaptor=$registry->get_adaptor( 'Multi', 'compara', 'SpeciesSet' );
  my $ref_all = $method_link_species_set_adaptor->fetch_all_by_method_link_type( $method_link_type );
  foreach my $mlssa (@$ref_all){
    print "name: ".$mlssa->name." method link type: ".$mlssa->method->type();
#    print "name: ".$mlssa->name." method link type: ".$mlssa->method_link_type; ## pre 70
#    print "" species set id:".$mlssa->species_set_id." source:".$mlssa->source."\n";
    my $specset=$species_set_adaptor->fetch_by_dbID( $mlssa->species_set_obj->dbID() );
#    my $specset=$species_set_adaptor->fetch_by_dbID( $mlssa->species_set_id );  ## pre 70
    print " species set:".$specset->name()."\n\n";
  }
}



