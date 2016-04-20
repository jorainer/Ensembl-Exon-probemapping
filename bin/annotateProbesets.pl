# This scripts retrieves annotation information of probe sets from ENSEMBL
# It can be used for different type of micro arrays from different vendors
my $script_version="0.0.1";

# loading additional modules
use DBI;
use strict;
use Getopt::Std;
use Bio::EnsEMBL::ApiVersion;
use Bio::EnsEMBL::Registry;

#### default settings:
my $array_name="HG-U133_Plus_2";
my $array_vendor="AFFY";
my $species="Human";
## connect to the public Ensembl database... is even faster than the local one.
my $user="anonymous";
my $host="ensembldb.ensembl.org";
my $mysql_port="5306";
my $pass="";

# Read option keys
my %option=();
getopts("c:h:s:v:",\%option);
if( $option{ h } ){
    print( "\nannotateProbesets.pl ".$script_version."\n" );
    print( "Retrieve all probe sets for a specified microarray along with their annotation to Ensembl genes.\n\n" );
    print( "usage: annotateProbesets.pl c:hs:v:\n" );
    print( "parameters:\n" );
    print( "-c (optional): the microarray/chip type. Defaults to HG-U133_Plus_2.\n" );
    print( "-h print this help.\n\n" );
    print( "-s (optional): the species. Defaults to human.\n" );
    print( "-v (optional): the microarray vendor. Defaults to AFFY.\n" )
    exit;
}
if( defined{ $option{ c } } ){
  $array_name=$option{ c };
}
if( defined{ $option{ s } } ){
  $species=$option{ s };
}
if( defined{ $option{ v } } ){
  $array_vendor=$option{ v };
}

## connecting to the Ensembl database.
my $registry='Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(-host => $host, -user => $user, -pass => $pass, -port => $mysql_port, -verbose => "1" );
my $ensembl_version=$registry->software_version;

my $info_string="# annotateProbesets.pl version $script_version querying Ensembl version ".$ensembl_version."; ".localtime()."\n# Array name: $array_name.\n# Array vendor: $array_vendor.\n# Species: $species\n";
my @probesets;

## if input file is defined we read the probe sets from there...
if( defined( $option{ i } ) ){
    my $infile=$option{ i };
    print "Reading probe sets from file: ".$infile." as source file !!!\n";
# load probe sets for annotation
    open(IN, "<$option{ i }") or die "cannot open file ".$option{ i }." for reading!\n";
    my $incount=0;
    while(<IN>){
	chomp;
	my $current_line = $_;
	$probesets[$incount] = $current_line;
	$incount++;
    }
    close( IN );
    $info_string=$info_string."# probe set ids read from file: ".$infile."\n";
}

if( defined( $option{ s } ) ){
    $species = $option{s};
}
my $outfile=$array_name."_Ensembl_".$ensembl_version.".txt";
print $info_string;
open( OUT, "> $outfile" ) or die "can't open file $outfile for writing!\n";
print OUT $info_string;
print OUT "probeset_id\tensembl_gene_id\tsymbol\tentrezid\texternal_name\tgene_biotype\tchromosome_name\tchromosome_strand\n";

## otherwise get all probe sets from the database...
my $array_adaptor = $registry->get_adaptor( $species, "funcgen", "array" );
my $array = $array_adaptor->fetch_by_name_vendor($array_name, $array_vendor);
my $tx_adaptor = $registry->get_adaptor( $species,"core","Transcript");
my $gene_adaptor = $registry->get_adaptor( $species,"core","Gene");
## fetch all probe sets.
@probesets = @{$array->get_all_ProbeSets()};
my $counter=0;
print "# querying ".scalar(@probesets)." probe sets.\n";
foreach my $current_probeset (@probesets){
    my @dbeList = @{$current_probeset->get_all_DBEntries()};
    my $transcript_string="";
    my %transcript_hash;
    my %gene_hash;
    my %symbol_hash;
    my %entrezid_hash;
    my %external_name_hash;
    my %chromosome_hash;
    my %chromosome_strand_hash;
    my %gene_biotype_hash;
    my $symbols_string="";
    my $entrezids_string="";
    my $genes_string="";
    my $external_names_string="";
    my $chromosomes_string="";
    my $chromosome_strands_string="";
    my $gene_biotypes_string="";
    if( scalar @dbeList > 0 ){
	foreach my $dbe (@dbeList){
	    ## determine which of the entries are Transcripts...
	    my $dbe_dbname = $dbe->dbname();
	    if ($dbe_dbname =~ /core_Transcript$/){
		my $gene = $gene_adaptor->fetch_by_transcript_stable_id($dbe->primary_id());
		#print " gene:".$gene->stable_id."";
		$gene_hash{ $gene->stable_id } = 0;
		## get symbol, entrezid, chromloc, description.
		#my $all_entries = $gene->get_all_DBLinks( "EntrezGene" );
		#foreach my $dbentry ( @{$all_entries} ){
		#    $symbol_hash{ $dbentry->display_id } = 0;
		#    $entrezid_hash{ $dbentry->primary_id } = 0;
		#}
		my $all_entries = $gene->get_all_DBEntries( "EntrezGene" );
		foreach my $dbentry ( @{$all_entries} ){
		    #print $dbentry->status."\n";
		    $symbol_hash{ $dbentry->display_id } = 0;
		    $entrezid_hash{ $dbentry->primary_id } = 0;
		}
		$external_name_hash{ $gene->external_name } = 0;
		$gene->transform("chromosome");
		$chromosome_hash{ $gene->slice->seq_region_name } = 0;
		$chromosome_strand_hash{ $gene->strand } = 0;
		$gene_biotype_hash{ $gene->biotype } = 0;
		$external_name_hash{ $gene->external_name } = 0;
	    }
	}
	$genes_string = join(";", keys %gene_hash);
	$symbols_string = join( ";", keys %symbol_hash );
	$entrezids_string = join( ";", keys %entrezid_hash );
	$chromosomes_string = join(";", keys %chromosome_hash);
	$chromosome_strands_string = join(";", keys %chromosome_strand_hash);
	$gene_biotypes_string = join(";", keys %gene_biotype_hash);
	$external_names_string = join( ";", keys %external_name_hash );
    }else{
#	print "no genes\n";
    }
    ## just some feedback...
    $counter++;
    if( ( $counter % 1000 ) == 0 ){
	print localtime()."> processed ".$counter." probe sets\n";
    }
    print OUT $current_probeset->name."\t".$genes_string."\t".$symbols_string."\t".$entrezids_string."\t".$external_names_string."\t".$gene_biotypes_string."\t".$chromosomes_string."\t".$chromosome_strands_string."\n";
    print ".".$counter;
}

close( OUT );
$registry -> disconnect_all();

