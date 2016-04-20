#!/usr/bin/perl
#####################################
my $script_version = "0.0.1";
##
# script to get probe sets matching/targeting genes defined in Ensembl.

#####################################
use IO::File;
use DBI;
use Getopt::Std;
use strict;
use warnings;
use Bio::EnsEMBL::ApiVersion;
use Bio::EnsEMBL::Registry;

my %option=();
getopts("a:c:hm:p:s:",\%option);
if( $option{ h } ){
    print( "\nannotate_chip.pl ".$script_version."\n" );
    print( "Retrieve all probe sets for all transcripts/genes defined in Ensembl.\n\n" );
    print( "usage: annotate_chip.pl a:c:hp:s:\n" );
    print( "parameters:\n" );
    print( "-a (optional): the minimum perfect alignment of a probe (defaults to 24, i.e. 25 nucleotides have to be aligned).\n" );
    print( "-c (optional): the microarray/chip type. Defaults to HG-U133_Plus_2.\n" );
    print( "-h print this help.\n\n" );
    print( "-m (optional): the maximal number of allowed mismatches for a probe. Defaults to 0." );
    print( "-p (optional): the minimum number of probes that have to be aligned within the exons of a transcript. Defaults to 9." );
    print( "-s (optional): the species. Defaults to human." );
    exit;
}

my $species="human";
my $chip_type="HG-U133_Plus_2";
my $min_probes=9;
my $min_probe_alignment=24;
my $max_mm=0;
my $probe_length=25;

if( defined( $option{ a } ) ){
  $min_probe_alignment=$option{ a };
}
if( defined( $option{ c } ) ){
  $chip_type = $option{ c };
}
if( defined( $option{ p } ) ){
  $min_probes = $option{ p };
}
if( defined( $option{ s } ) ){
  $species = $option{ s };
}

## Ensembl settings.
my $ensembl_user="anonymous";
my $ensembl_host="ensembldb.ensembl.org";
my $ensembl_password="";
my $ensembl_port="5306";


## check Ensembl version...
my $api_version="".software_version()."";


my $out_file = "".$chip_type."-Ensembl-".$api_version.".txt";
my $infostring = "# annotate_chip.pl version $script_version\n# species: $species\n# chip_type: $chip_type\n# min_probes: $min_probes\n# min_probe_alignment: $min_probe_alignment\n# max mismatches: $max_mm\n# Ensembl version: $api_version\n";

print( $infostring );

my $registry = 'Bio::EnsEMBL::Registry';
$registry->set_reconnect_when_lost();
$registry->load_registry_from_db(-host => $ensembl_host, -user => $ensembl_user, -pass => $ensembl_password, -verbose => "1" );
my $gene_adaptor = $registry->get_adaptor( $species, "core", "gene" );
my $probe_feature_adaptor = $registry->get_adaptor( $species, "funcgen", "ProbeFeature" );
my $aa = $registry->get_adaptor( $species, "funcgen", "Array" );
my $array = $aa->fetch_by_name_vendor( $chip_type, 'AFFY');
open( OUT , ">$out_file") or die "cannot open file ".$out_file." for writing!\n";
print OUT $infostring;
print OUT "probeset_id\tno_probes\ttranscript_id\tgene_id\tgene_name\tsymbol\tentrezgene\tgene_biotype\tchromosome_name\tcoord_system\tstrand\n";
## get all gene ids defined in the database:
## 1) get gene,
## 2) get all transcripts
### 3) get probe sets.
print "start fetch...";
my @gene_ids = @{$gene_adaptor->list_stable_ids()};
print "done\n";
my $counta = 0;
foreach my $gene_id ( @gene_ids ){
  $counta++;
  print "$gene_id\n";
  if( ($counta % 100) == 0 ){
    print localtime()."> Processed $counta of ".scalar @gene_ids." genes.\n";
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
    if( !$gene->is_known | !$gene->is_current ){
      next;
    }
    my $coord_system = $gene->coord_system_name;
    my $chromosome_name = $gene->slice->seq_region_name;
    my $strand = $gene->strand;
    my $gene_external_name= $gene->external_name;
    my $gene_biotype = $gene->biotype;
    my $symbol="";
    my $entrezgene="";
    #my $all_entries = $gene->get_all_DBLinks( "EntrezGene" );
    #my %symbol_hash=();
    #my %entrezgene_hash=();
    #foreach my $dbe ( @{$all_entries} ){
    #  $symbol_hash{ $dbe->display_id } = "something";
    #  $entrezgene_hash{ $dbe->primary_id } = "something";
    #}
    #$symbol = join( ";", sort keys %symbol_hash );
    #$entrezgene = join( ";", sort keys %entrezgene_hash );
    my @transcripts = @{ $gene->get_all_Transcripts };
    ## looping through the transcripts
    foreach my $transcript ( @transcripts ){
#      if( $do_transform==1 ){
#	$transcript = $transcript->transform( "chromosome" );  ## only transform if gene can be mapped to a valid chromosome.
#      }
      my $transcript_id=$transcript->stable_id;
      ## OK, we're doing this terribly slow. we're acutally going through each exon and fetching the data from there.
      my %transcript_probesets=();	## putting the probeset ids and counts per probeset id into this hash.
      my @exons = @{ $transcript->get_all_Exons() };
      foreach my $exon ( @exons ){
	if( $do_transform==1 ){
	  $exon = $exon->transform( "chromosome" );
	}
	my $slice = $exon->feature_Slice;
	if( defined $slice ){
	  ##my @features=();
	  my @features = @{ $probe_feature_adaptor->fetch_all_by_Slice_Array( $slice, $array )};
	  ##my @features = @{ $probe_feature_adaptor->fetch_all_by_Slice_array_vendor( $slice, $chip_type, "AFFY" )};
	  if( scalar( @features ) > 0 ){
	    foreach my $feature (@features){
	      ## ok if alignment is 24 and no mismatches add it to a hash...
	      if( $feature->mismatchcount <= $max_mm && ( $feature->end - $feature->start ) >= $min_probe_alignment && ( $feature->end - $feature->start ) < $probe_length && $feature->strand == $strand){
		my $ps = $feature->probe->probeset->name;
		if( exists( $transcript_probesets{ $ps } ) ){
		  my $dummy = $transcript_probesets{ $ps };
		  $dummy++;
		  $transcript_probesets{ $ps } = $dummy;
		}
		else{
		  $transcript_probesets{ $ps } = 1;
		}
	      }
	    }
	  }
	}
      }
      ## writing the output... if we do have any probe sets left...
      if( scalar( keys %transcript_probesets ) > 0 ){
	print "got ".scalar( keys %transcript_probesets )." probe sets for transcript ".$transcript_id.":";
	foreach my $key ( keys %transcript_probesets ){
	  if( $transcript_probesets{ $key } >= $min_probes ){
	    print " gotcha! ".$key.";";
	    print OUT $key."\t".$transcript_probesets{ $key }."\t".$transcript_id."\t".$gene_id."\t".$gene_external_name."\t".$symbol."\t".$entrezgene."\t".$gene_biotype."\t".$chromosome_name."\t".$coord_system."\t".$strand."\n";
	  }else{
	    print " ".$key." only ".$transcript_probesets{ $key }." probes;"
	  }
	}
	print "\n";
      } ## end if
      else{
      }
    }   ## end foreach my $transcript
  }     ## end if( defined $orig_gene )
}

close( OUT );


#####
## here we read the above created file and re-format it that each row contains a
## single probe set with all other elements eventually pasted and separated by ;
sub generate_annotation_by_probeset{
}


## other version... get by probe set.
#sub annotate_via_probeset_id{
#  my $probeset_adaptor = $registry->get_adaptor( $species, 'funcgen', 'ProbeSet' );
#  my @probesets = @{ $probeset_adaptor -> fetch_all_by_Array( $array ) };
#  foreach my $probeset (@probesets){
#    $counta++;
#    if( ($counta % 100) == 0 ){
#      print localtime()."> Processed $counta of ".scalar @probesets." probe sets.\n";
#    }
#  }
#}

