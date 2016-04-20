package CustomCDF::Transcript;
use strict;
use warnings;
use Carp;
use CustomCDF::Feature;
our $VERSION = "0.3.0";

our @ISA = qw( CustomCDF::Feature );

=head1 NAME

Transcript - A transcript representation

=head1 SYNOPSIS

This is a representation of a transcript. Transcripts can be retrieved from an alignment and annotation database using:

my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost",
						  username=>"anonuser",
						  password=>"",
						  dbname=>"homo_sapiens_hugene_cdna_74"
				     );
my $gene_adaptor = $db_adaptor->get_gene_adaptor();
my $transcript = $gene_adaptor->fetch_transcript( transcript_id=>"ENST00000373020", load_exons=>1 );
print "call print on transcript:\n";
$transcript->print( header=>1 );

=head1 DESCRIPTION

Is a representation of a transcript.
Attributes are:

=item transcript_id                    : the (Ensembl) transcript id

=item start/transcript_chrom_start     : the genomic start position of the transcript

=item end/transcript_chrom_end         : the genomic end position.

=item strand/transcript_chrom_strand   : the strand on which the transcript is encoded.

=item transcript_coding_chrom_start    : the genomic position of the start of the transcript's coding region.

=item transcript_coding_chrom_end      : the genomic position of the end of the transcript's coding region.

=item transcript_biotype               : the transcript's biotype.

=item chromosome_name                  : the name of the chromosome on which the transcript is encoded.

=item coord_system                     : the coordinate system of the gene; for genes encoded on.

=item exons                            : an array with exons

=item adaptor                          : the CustomCDF:GeneAdaptor.

CustomCDF:Transcript extends CustomCDF::Feature

=head2 Methods

=cut

sub new{
  my($class, %args) = @_;
  my $self = bless({}, $class);
  $self->_init(%args);
  return $self;
}

sub _init{
  my ($self, %args)=@_;
  my $transcript_id = exists $args{transcript_id} ?  $args{transcript_id} : "NA";
  my $transcript_coding_chrom_start = exists $args{transcript_coding_chrom_start} ?  $args{transcript_coding_chrom_start} : -1;
  my $transcript_coding_chrom_end = exists $args{transcript_coding_chrom_end} ?  $args{transcript_coding_chrom_end} : -1;
  my $transcript_biotype = exists $args{transcript_biotype} ?  $args{transcript_biotype} : "NA";
  #my $chromosome_name = exists $args{chromosome_name} ?  $args{chromosome_name} : "NA";
  #my $coord_system = exists $args{coord_system} ?  $args{coord_system} : "NA";
  my $exons = [];
  if( exists $args{exons} ){
    $exons=$args{exons};
    $self->exons( $exons );
  }
  $self->transcript_id( $transcript_id );
  $self->transcript_coding_chrom_start( $transcript_coding_chrom_start );
  $self->transcript_coding_chrom_end( $transcript_coding_chrom_end );
  $self->transcript_biotype( $transcript_biotype );
  ## things passed to Feature:
  # start, end, strand, coord_system, chromosome_name, adaptor
  my $transcript_chrom_start = exists $args{transcript_chrom_start} ?  $args{transcript_chrom_start} : -1;
  if( !exists $args{ start } ){
    $args{ start } = $transcript_chrom_start;
  }
  my $transcript_chrom_end = exists $args{transcript_chrom_end} ?  $args{transcript_chrom_end} : -1;
  if( !exists $args{ end } ){
    $args{ end } = $transcript_chrom_end;
  }
  my $transcript_chrom_strand = exists $args{transcript_chrom_strand} ?  $args{transcript_chrom_strand} : 0;
  if( !exists $args{ strand } ){
    $args{ strand } = $transcript_chrom_strand;
  }
  $self->SUPER::_init(%args);
}

=head3 transcript_id

Getter/setter for the transcript_id attribute.

=cut
sub transcript_id{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{transcript_id} = $x;
  }
  return $self->{transcript_id};
}


=head3 transcript_coding_chrom_start

Getter/setter for the (chromosomal) start position of the transcript's coding region.

Note: returns -1 for non-coding transcripts!

=cut
sub transcript_coding_chrom_start{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    if( !defined $x ){
      $x = -1;
    }
    $self->{transcript_coding_chrom_start} = $x;
  }
  return $self->{transcript_coding_chrom_start};
}

=head3 transcript_coding_chrom_end

Getter/setter for the (chromosomal) end position of the transcript's coding region.

Note: returns -1 for non-coding transcripts!

=cut
sub transcript_coding_chrom_end{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    if( !defined $x ){
      $x = -1;
    }
    $self->{transcript_coding_chrom_end} = $x;
  }
  return $self->{transcript_coding_chrom_end};
}


=head3 transcript_biotype

Getter/setter for the biotype of the transcript.

=cut
sub transcript_biotype{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{transcript_biotype} = $x;
  }
  return $self->{transcript_biotype};
}


=head3 exons

Getter/setter for the exons parameter. Exons will be empty by default, unless the transcript is fetchted with the additional parameter "load_exons=>1". Alternatively, the exons can be loaded later using the "load_exons" method.

=cut
sub exons{
  my $self = shift;
  if( @_ ) {
    my $exons = shift;
    croak "Error: illegal value for exons; should be a reference to an array!" unless ref $exons eq "ARRAY";
    my @exons_array = @{ $exons };
    foreach my $the_exon ( @exons_array ){
      if( !eval { $the_exon->isa( 'CustomCDF::Exon' ) } ){
	croak "Error: exons should be a reference to CustomCDF::Exon objects!";
      }
    }
    $self->{exons}=$exons;
  }
  return @{$self->{exons}||[]};
}


=head3 load_exons

Retrieves (and returns) all exons for the transcript from the database and stores them into the "exons" attribute of the transcript.

=cut
sub load_exons{
  my $self = shift;
  my @exons =$self->fetch_exons();
  $self->exons( \@exons );
  return @exons;
}


=head3 fetch_exons

Retrieves all exons for the transcript. Returns an array of CustomCDF::Exon, exons are ordered by their chromosomal start position. Note that this call will always result in a database query. Alternatively, exons can be loaded using the load_exons method, which will store them locally in the exons parameter of the CustomCDF::Transcript.

=cut
sub fetch_exons{
  my $self = shift;
  return ($self->adaptor()->fetch_exons_for_transcript( transcript_id=>$self->transcript_id ) );
}

=head3 fetch_gene

Retrieves the Gene object of the transcript's gene.

=cut
sub fetch_gene{
  my $self=shift;
  return( $self->adaptor()->fetch_gene_by_transcript( transcript_id=>$self->transcript_id ) );
}


=head3 to_string

Pastes all values of the object to a tab-delimited string.

=cut
sub to_string{
  my $self=shift;
  my $have_exons = scalar $self->exons();
  my $the_string=$self->transcript_id."\t".$self->transcript_biotype."\t".$self->start."\t".$self->end."\t".$self->transcript_coding_chrom_start."\t".$self->transcript_coding_chrom_end."\t".$self->chromosome_name."\t".$self->strand."\t".$self->coord_system;
  if( $have_exons > 0 ){
    $the_string = $the_string."\t".$have_exons;
  }
  return( $the_string );
}

=head3 print

Prints the values for the probe.
If an argument "header" is provided it will also print a header.

=cut
sub print{
  my $self=shift;
  my %args = @_;
  if( exists( $args{ header } ) ){
    my $header = "transcript_id\ttranscript_biotype\tstart\tend\ttranscript_coding_chrom_start\ttranscript_coding_chrom_end\tchromosome_name\tstrand\tcoord_system";
    my $have_exons = scalar $self->exons();
    if( $have_exons > 0 ){
      $header = $header."\tn_exons";
    }
    print $header."\n";
  }
  ## just printing the data... tab separated...
  print $self->to_string."\n";
}

=head3 length

Returns the length of the transcript in nucleotides.

=cut
sub length{
  my $self=shift;
  my $length=0;
  my @exons = $self->exons();
  foreach my $exon (@exons){
    $length=$length+$exon->length;
  }
  return $length;
}



1;

