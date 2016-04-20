## test script for ProbeAdaptor, Probe, Alignment.

use CustomCDF::ProbeAdaptor;
use CustomCDF::DBAdaptor;
use CustomCDF::Probe;

# my $probe_adaptor = CustomCDF::ProbeAdaptor->new( host=>"localhost",
# 						  username=>"anonuser",
# 						  password=>"",
# 						  dbname=>"homo_sapiens_hugene_cdna_74"
# 				     );
my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost",
						  username=>"anonuser",
						  password=>"",
						  dbname=>"homo_sapiens_hugene_cdna_74"
				     );
my $probe_adaptor = $db_adaptor->get_probe_adaptor();

# my $dbcon=$probe_adaptor->dbcon();
# my $query = $dbcon->prepare( "select * from probe_to_exon_transcript where chromosome_strand=1 limit 1" );
# $query->execute();
# my %hash = %{$query->fetchrow_hashref};
# print "the chromosome_strand: ".$hash{chromosome_strand}."\n";
# if( $hash{ chromosome_strand }==1 ){
#   print "chromosome_strand==1\n";
# }else{
#   print "chromosome_strand!=1\n";
# }

# my $dbcon=$probe_adaptor->dbcon();
# my $query = $dbcon->prepare( "select * from probe_alignments where seq_strand='-' limit 1" );
# $query->execute();
# my %hash = %{$query->fetchrow_hashref};
# print "the chromosome_strand: ".$hash{seq_strand}."\n";
# if( $hash{ seq_strand } eq "+" ){
#   print "chromosome_strand eq +\n";
# }else{
#   print "chromosome_strand ne +\n";
# }



#print "fetching all probes\n";
#my @probes = $probe_adaptor->fetch_all_probes();
#print "got: ".scalar(@probes)." probes\n";


print "\nfetching probes which map to a single gene and do not have more than 1 alignment to the genome (also considering a missmatch)\n";
my @probes = $probe_adaptor->fetch_all_probes( nr_gene_map=>"=1", nr_gene_map_mm=>"<=0", nr_chrom_map=>"<=1", nr_chrom_map_mmall=>"=0" );
print "got: ".scalar(@probes)." probes\n";
print "some informations on the first probe:\n";
print "probe sequence: ".$probes[0]->sequence."\n";
print "does the probe have an adaptor?\n";
print "answer:".$probes[0]->have_adaptor()."\n";
print "printing the entries of the first probe:\n";
$probes[ 0 ]->print( header=>1 );

## that's perfectly working. Stores a reference rather than a copy!
# ### try to close the dbconnection from one and then test it again...
# print "closing the dbcon on adaptor of first probe.";
# $probes[0]->adaptor->dbcon->disconnect();
# print "and now get the dbcon from the next probe. If a copy was stored we should be able to execute stuff on that...\n";
# my $con <- $probes[1]->adaptor->dbcon();
# my $query = $con->prepare( "show tables" ) or die $con->errstr();
# $query->execute or die $con->errstr();
# while (my @row = $query->fetchrow_array){
#   print "got ".$row[0]."\n";
# }

## ok, now we are fetching probe alignments for the first probe.
print "Fetching all alignments to cdna for the first probe\n";
$probes[ 0 ]->load_alignments( seq_type=>"='cdna'" );
$probes[ 0 ]->print( header=>1 );
print "Printing all alignments:\n\n";
my @aligns = $probes[ 0 ]->alignments;
my $do_header=1;
foreach my $align ( @aligns ){
  if( $do_header==1 ){
    $align->print( header=>1 );
    $do_header=0;
  }else{
    $align->print();
  }
}

## fetch genomic alignments for probe.
print "\nFetching all genomic alignments (without mismatch) for the probe (513856) targeting a gene on + strand (with splice junction alignment).\n";
@aligns = $probe_adaptor->fetch_genomic_alignments_for_probe( probe_id=>"513856" );
my $do_header=1;
foreach my $align ( @aligns ){
  if( $do_header==1 ){
    $align->print( header=>1 );
    $do_header=0;
  }else{
    $align->print();
  }
}



print "\nFetching all genomic alignments (without mismatch) for the probe (220661) targeting a gene on the - strand.\n";
@aligns = $probe_adaptor->fetch_genomic_alignments_for_probe( probe_id=>"220661" );
my $do_header=1;
foreach my $align ( @aligns ){
  if( $do_header==1 ){
    $align->print( header=>1 );
    $do_header=0;
  }else{
    $align->print();
  }
}


print "\nFetching a probe from the database and pre-loading its alignments against cdna.\n";
my $probe = $probe_adaptor->fetch_probe( probe_id=>"156902", load_alignments=>1, seq_type=>"='cdna'" );
$probe->print( header=>1, alignment=>1 );

print "\nLoading alignments against cdna and/or chromosomes for a probe from the database.\n";
$probe->load_alignments( );
$probe->print( header=>1, alignment=>1 );

print "\nFetching a probe from the database and pre-loading its genomic alignments.\n";
my $probe = $probe_adaptor->fetch_probe( probe_id=>"156902", load_genomic_alignments=>1 );
$probe->print( header=>1, alignment=>1 );



##
print "\n\n------------------------------------------\n Exon related stuff.\n";
print "fetch all probes for a given exon:\n";
@probes = $probe_adaptor->fetch_probes_for_exon( exon_id=>"ENSE00003478088", all=>1 );
my $do_header=1;
foreach my $probe ( @probes ){
  if( $do_header==1 ){
    $probe->print( header=>1 );
    $do_header=0;
  }else{
    $probe->print();
  }
}
print "\nfetch good probes with alignment for a given exon (ENSE00003478088):\n";
@probes = $probe_adaptor->fetch_probes_for_exon( exon_id=>"ENSE00003478088", load_alignments=>1 );
my $do_header=1;
foreach my $probe ( @probes ){
  if( $do_header==1 ){
    $probe->print( header=>1, alignment=>1 );
    $do_header=0;
  }else{
    $probe->print( alignment=>1 );
  }
}

print "\nfetch good probes with genomic alignments for a given exon (ENSE00003478088):\n";
@probes = $probe_adaptor->fetch_probes_for_exon( exon_id=>"ENSE00003478088", load_genomic_alignments=>1 );
my $do_header=1;
foreach my $probe ( @probes ){
  if( $do_header==1 ){
    $probe->print( header=>1, alignment=>1 );
    $do_header=0;
  }else{
    $probe->print( alignment=>1 );
  }
}


##
print "\n------------------------------------------\n Transcript related stuff.\n";
print "\nFetch all probes for a given transcript (ENST00000311549). First without probe_alignments.\n";
@probes = $probe_adaptor->fetch_probes_for_transcript( transcript_id=>"ENST00000311549", all=>1 );
my $do_header=1;
foreach my $probe ( @probes ){
  if( $do_header==1 ){
    $probe->print( header=>1, alignment=>1 );
    $do_header=0;
  }else{
    $probe->print( alignment=>1 );
  }
}

print "\nFetch all probes for a given transcript (ENST00000311549). Load also all alignments.\n";
@probes = $probe_adaptor->fetch_probes_for_transcript( transcript_id=>"ENST00000311549", all=>1, load_alignments=>1 );
my $do_header=1;
foreach my $probe ( @probes ){
  if( $do_header==1 ){
    $probe->print( header=>1 );
    $do_header=0;
  }else{
    $probe->print( );
  }
}

print "\nFetch good probes for a given transcript (ENST00000395766) with all their genomic alignments (note: we do have junction probes for this one).\n";
@probes = $probe_adaptor->fetch_probes_for_transcript( transcript_id=>"ENST00000395766", load_genomic_alignments=>1 );
my $do_header=1;
foreach my $probe ( @probes ){
  if( $do_header==1 ){
    $probe->print( header=>1, alignment=>1 );
    $do_header=0;
  }else{
    $probe->print( alignment=>1 );
  }
}


print "\nSome more tests on transcripts and their probes.\n";
@probes = $probe_adaptor->fetch_probes_for_transcript( transcript_id=>"ENST00000159647" );
my $do_header=1;
foreach my $probe ( @probes ){
  if( $do_header==1 ){
    $probe->print( header=>1, alignment=>1 );
    $do_header=0;
  }else{
    $probe->print( alignment=>1 );
  }
}

## do we have any gapped alignments at al for "good" probes?
# my $dbcon=$probe_adaptor->dbcon();
# my $query = $dbcon->prepare( "select distinct transcript_id from probe_to_exon_transcript where is_junction_alignment=1" );
# $query->execute();
# while ( my $hashref = $query->fetchrow_hashref ){
#   my %rowhash = %{$hashref};
#   my $transcript_id = $rowhash{ transcript_id };
#   print "checking transcript ".$transcript_id.".";
#   my @probes = $probe_adaptor->fetch_probes_for_transcript( transcript_id=>$transcript_id, load_genomic_alignments=>1 );
#   foreach my $probe ( @probes ){
#     print ".";
#     ## check whether the probe has a gapped alignment...
#     my @algs = $probe->alignments();
#     foreach my $alg ( @algs ){
#       if( $alg->is_gapped() == 1 ){
# 	print "\n hooray! got a exon junction probe!\n";
# 	$probe->print( alignment=>1 );
#       }
#     }
#   }
#   print "done\n";
# }


print "\n------------------------------------------\n Gene related stuff.\n";
print "\nFetch all probes for a given gene (ENSG00000122026). First without probe_alignments.\n";
@probes = $probe_adaptor->fetch_probes_for_gene( gene_id=>"ENSG00000122026", all=>1 );
my $do_header=1;
foreach my $probe ( @probes ){
  if( $do_header==1 ){
    $probe->print( header=>1, alignment=>1 );
    $do_header=0;
  }else{
    $probe->print( alignment=>1 );
  }
}


print "\nFetch all probes for a given gene (ENSG00000122026) with genomic alignments.\n";
@probes = $probe_adaptor->fetch_probes_for_gene( gene_id=>"ENSG00000122026", load_genomic_alignments=>1, all=>1 );
my $do_header=1;
foreach my $probe ( @probes ){
  if( $do_header==1 ){
    $probe->print( header=>1 );
    $do_header=0;
  }else{
    $probe->print( );
  }
}


print "\nFetch good probes for a given gene (ENSG00000122026) with alignments to chromosome.\n";
@probes = $probe_adaptor->fetch_probes_for_gene( gene_id=>"ENSG00000122026", load_alignments=>1, seq_type=>"='chromosome'" );
my $do_header=1;
foreach my $probe ( @probes ){
  if( $do_header==1 ){
    $probe->print( header=>1, alignment=>1 );
    $do_header=0;
  }else{
    $probe->print( alignment=>1 );
  }
}


## just a timing comparison
# print "\n\nTiming: probes for exon with genomic alignments.\n";
# print ">".localtime()."\n";
# @probes = $probe_adaptor->fetch_probes_for_gene( gene_id=>"ENSG00000122026", load_genomic_alignments=>1 );
# print ">".localtime()."\n";


# print "\n\nTiming: probes for exon with alignments against chromosome.\n";
# print ">".localtime()."\n";
# @probes = $probe_adaptor->fetch_probes_for_gene( gene_id=>"ENSG00000122026", load_alignments=>1, seq_type=>"='chromosome'" );
# print ">".localtime()."\n";

print "\n\nlooks fine\n";
