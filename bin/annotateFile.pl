#######
## 14.02.2011
my $script_version_="0.0.3";
##### description
## just a very cool perl script to annotate a tabulatir delimited text file with some additional columns from an
## ensembl database.
# version 0.0.2: 02.02.2011: add chromosome name, strand and gene biotype
# version 0.0.3: 14.02.2011: use join instead of "concatenate"
# paramters:
# s: species, defaults to human.
# g: column in the input file containing the ensembl gene IDs, defaults to gene_id

#use lib '/home500/ensembl/API/current_API/ensembl/ensembl/modules';
#use lib '/home500/ensembl/temp/52_36n/ensembl/ensembl/modules';
use strict;
use Getopt::Std;
use Bio::EnsEMBL::Registry;

my %option=();
my $ensembl_database="core";
my $species="human";
my $gene_id_column="gene_id";
my $input_file;
my $output_file;
getopts("s:g:i:o:h",\%option);
if( defined( $option{ h } ) ){
  print( "\nannotateFile.pl version ".$script_version_.".\n" );
  print( "Retrieve various annotations for Ensembl gene IDs provided in the input file using the Ensembl Perl API. Multiple values for a gene will be pasted into a single value, separated by a semicolon (;).\n\n");
  print( "Usage: perl annotateFile.pl -i:o:[g:hs:]\n\n" );
  print( "Parameters:\n" );
  print( "-i (required): the input file, i.e. a tabulator delimited text file that contains the Ensembl gene IDs, one ID per row.\n" );
  print( "-o (required): the output file to which the annotation should be written. The output file will contain all columns from the input file with the annotations in additional columns.\n" );
  print( "-g (optional): the name of the column in the input file containing the Ensembl gene IDs (defaults to gene_id).\n" );
  print( "-h : print this help and exit\n" );
  print( "-s (optional): the species; defaults to human.\n" );
  exit 0;
}
if( defined( $option{ s } ) ){
  $species=$option{ s };
}
if( defined( $option{ g } ) ){
  $gene_id_column=$option{ g };
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

my $user="anonuser";
my $host="madb.i-med.ac.at";
#my $user="anonymous";
#my $host="ensembldb.ensembl.org";
#my $port="5306";
my $port="3306";
my $pass="";
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(-host => $host, -user => $user, -pass => $pass, -port => $port, -verbose => "1" );

my $gene_adaptor = $registry->get_adaptor( $species, $ensembl_database, "gene" );


# Define a helper subroutine to print DBEntries
sub print_DBEntries
{
    my $db_entries = shift;

    foreach my $dbe ( @{$db_entries} ) {
        printf "\tXREF %s (%s) primary id %s\n", $dbe->display_id(), $dbe->dbname(), $dbe->primary_id;
    }
}

#sub concatenate{
#  my @the_keys = shift;
#  my $string="";
#  my $counta=0;
#  foreach my $key ( @the_keys ){
#    if( $counta==0 ){
#      $string=$key;
#      $counta = 1;
#    }
#    else{
#      $string=$string.";".$key;
#    }
#  }
#  return $string;
#}

my $gene_id_idx = -1;
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
	if( $colname eq $gene_id_column ){
	  $gene_id_idx = $colcounter;
	}
	$colcounter=$colcounter+1;
      }
      if( $gene_id_idx < 0 ){
	die( "column $gene_id_column not found in header of input file!\n" );
      }
      $header=1;
      print OUT $_."\tsymbol\tentrezgene\tgene_name\tgene_biotype\tchromosome_name\tchromosome_strand\n";
#      print OUT $_."\tgene_biotype\tchromosome_name\tchromosome_strand\n";
    }
    else{
      my %symbol_hash=();
      my %entrezgene_hash=();
#      print "gene_id: ".$line[ $gene_id_idx ]."\n";
      my $symbol="";
      my @symbols = ();
      my $entrezgene="";
      my @entrezgenes = ();
      my $external_name="";
      my @external_names = ();
      my $gene_biotype="";
      my @gene_biotypes = ();
      my $chromosome_name="";
      my @chromosome_names = ();
      my $chromosome_strand="";
      my @chromosome_strands = ();
      ## now we allow ; separated multiple gene ids...
      my @geneids = split( /;/, $line[ $gene_id_idx ] );
      foreach my $gene_id (@geneids){
	my $gene = $gene_adaptor->fetch_by_stable_id( $gene_id );
	if( defined( $gene ) ){
	  my $all_entries = $gene->get_all_DBLinks( "EntrezGene" );
	  foreach my $dbe ( @{$all_entries} ){
	    $symbol_hash{ $dbe->display_id } = "something";
	    $entrezgene_hash{ $dbe->primary_id } = "something";
	  }
	  push (@symbols, sort keys %symbol_hash );
	  push (@entrezgenes, sort keys %entrezgene_hash );
	  push (@external_names, $gene->external_name);
	  $gene->transform("chromosome");
	  push (@chromosome_names, $gene->slice->seq_region_name );
	  push (@gene_biotypes, $gene->biotype );
	  push (@chromosome_strands, $gene->strand );
	}
	else{
	  print "gene $gene_id not found in database\n";
	}
      }
      ## ok, and now collapsing all the values I've got... making unique if possible.
      $symbol = join( ";", makeUnique( @symbols ) );
      $entrezgene = join( ";", makeUnique( @entrezgenes ));
      $external_name = join( ";", makeUnique( @external_names ) );
      $gene_biotype = join( ";", makeUnique( @gene_biotypes ) );
      $chromosome_name = join( ";", makeUnique( @chromosome_names ) );
      $chromosome_strand = join( ";", makeUnique( @chromosome_strands ) );
      print OUT $_."\t".$symbol."\t".$entrezgene."\t".$external_name."\t".$gene_biotype."\t".$chromosome_name."\t".$chromosome_strand."\n";
#      print OUT $_."\t".$gene_biotype."\t".$chromosome_name."\t".$chromosome_strand."\n";
#      print "gene: $line[ $gene_id_idx ] symbol: $symbol entrezgene: $entrezgene\n";
    }
  }
}



close( IN );
close( OUT );



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

