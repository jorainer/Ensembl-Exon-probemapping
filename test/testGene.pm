## test script for GeneAdaptor, Gene, Transcript and Exon.

use CustomCDF::GeneAdaptor;
use CustomCDF::DBAdaptor;
use CustomCDF::Gene;

my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost",
						  username=>"anonuser",
						  password=>"",
						  dbname=>"homo_sapiens_hugene_cdna_74"
				     );
my $gene_adaptor = $db_adaptor->get_gene_adaptor();


## testing some gene stuff.
print "\n--------------------------------\n> testing gene stuff\n";
my $gene = $gene_adaptor->fetch_gene( gene_id=>"ENSG00000000003" );
print "call print on gene:\n";
$gene->print( header=>1 );
print "\nload also all transcripts for that gene:\n";
$gene = $gene_adaptor->fetch_gene( gene_id=>"ENSG00000000003", load_transcripts=>1 );
$gene->print( header=>1 );

print "\nprinting transcript information:\n";
my $do_header=1;
foreach my $transcript ( $gene->transcripts() ){
  if( $do_header==1 ){
    $do_header=0;
    $transcript->print( header=>1 );
  }else{
    $transcript->print();
  }
}

print "\nload all transcripts and exons for that gene:\n";
$gene = $gene_adaptor->fetch_gene( gene_id=>"ENSG00000000003", load_exons=>1 );
$gene->print( header=>1 );

print "\nprinting transcript information:\n";
my $do_header=1;
foreach my $transcript ( $gene->transcripts() ){
  if( $do_header==1 ){
    $do_header=0;
    $transcript->print( header=>1 );
  }else{
    $transcript->print();
  }
}

print "\nfetch all transcripts for that gene (without exons)\n";
$gene->load_transcripts();
my $do_header=1;
foreach my $transcript ( $gene->transcripts() ){
  if( $do_header==1 ){
    $do_header=0;
    $transcript->print( header=>1 );
  }else{
    $transcript->print();
  }
}



## loading all exons for a given gene.
print "\nload all exons for that gene:\n";
my @exons = $gene->fetch_exons();
$do_header = 1;
foreach my $exon ( @exons ){
  if( $do_header==1 ){
    $do_header=0;
    $exon->print( header=>1 );
  }else{
    $exon->print();
  }
}


## testing some transcript stuff.
print "\n--------------------------------\n> testing transcript stuff\n";
my $transcript = $gene_adaptor->fetch_transcript( transcript_id=>"ENST00000373020", load_exons=>1 );
print "call print on transcript:\n";
$transcript->print( header=>1 );
print "print the pre-loaded exons:\n";
@exons = $transcript->exons();
$do_header = 1;
foreach my $exon ( @exons ){
  if( $do_header==1 ){
    $do_header=0;
    $exon->print( header=>1 );
  }else{
    $exon->print();
  }
}
print "\nThe length of the transcript is:".$transcript->length."\n";


print "\nload all exons for that transcript:\nNote: here we are first loading the plain transcript without exons and are subsequently fetching the exons for that transcript.\n";
my $transcript = $gene_adaptor->fetch_transcript( transcript_id=>"ENST00000373020" );
@exons = $transcript->fetch_exons();
$do_header = 1;
foreach my $exon ( @exons ){
  if( $do_header==1 ){
    $do_header=0;
    $exon->print( header=>1 );
  }else{
    $exon->print();
  }
}





## testing some stuff on exons:
print "\n--------------------------------\n> testing exon stuff\n";
my $exon = $gene_adaptor->fetch_exon( exon_id=>"ENSE00000327880" );
print "call print on exon:\n";
$exon->print( header=>1 );

my $exon = $gene_adaptor->fetch_exon( exon_id=>"ENSE00001855382" );
print "call print on exon:\n";
$exon->print( header=>1 );








