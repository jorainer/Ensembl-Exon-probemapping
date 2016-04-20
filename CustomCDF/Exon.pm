package CustomCDF::Exon;
use strict;
use warnings;
use Carp;
use CustomCDF::Feature;
our $VERSION = "0.3.0";

our @ISA = qw( CustomCDF::Feature );

=head1 NAME

Exon - A exon representation

=head1 SYNOPSIS

This object represents an exon of a transcript. Exons can be fetched from an alignment and annotation database via:

my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost",
						  username=>"anonuser",
						  password=>"",
						  dbname=>"homo_sapiens_hugene_cdna_74"
				     );
my $gene_adaptor = $db_adaptor->get_gene_adaptor();
my $exon = $gene_adaptor->fetch_exon( exon_id=>"ENSE00000327880" );
print "call print on exon:\n";
$exon->print( header=>1 );

=head1 DESCRIPTION

Is a representation of an Exon.
Attributes are:

=item exon_id          : the (Ensembl) exon id

=item start : the chromosomal start position of the exon.

=item end   : the chromosomal end position of the exon.

=item chromosome_name  : the name of the chromosome on which the transcript is encoded.

=item strand           : the strand on which the exon is encoded.

=item coord_system     : the coordinate system of the gene.

=item transcripts      : array with transcript ids for that exon.

=item adaptor          : the CustomCDF::GeneAdaptor

CustomCDF::Exon extends CustomCDF::Feature

=head2 Methods

=cut
sub new{
  my $class = shift;
  #my($class, %args) = @_;
  my $self = bless({}, $class);
  $self->_init(@_);
  return $self;
}

sub _init{
  my ($self, %args)=@_;
  ## pick all from the $args that we need to store in the Feature, i.e. start, end, chromosome_name, strand and coord_system
  my $exon_id = exists $args{ exon_id } ? $args{ exon_id } : "NA";
  my $start = exists $args{ exon_chrom_start } ? $args{ exon_chrom_start } : -1;
  my $end = exists $args{ exon_chrom_end } ? $args{ exon_chrom_end } : -1;
  my $strand = exists $args{ transcript_chrom_strand } ? $args{ transcript_chrom_strand } : 0;
  if( exists $args{ strand } ){
    $strand = $args{ strand };
  }
  ## chromosome_name and coord_system are expected to be named like that...
  $args{ start } = $start;  ## eventually overwriting existing entries... that's ok.
  $args{ end } = $end;
  $args{ strand } = $strand;
  $self->exon_id( $exon_id );
  $self->SUPER::_init(%args);
  ## setting all exon related stuff.
}

=head3 exon_id

Getter/Setter for the exon_id.

=cut
sub exon_id{
  my $self = shift;
  if( @_ ){
    my $x = shift;
    $self->{exon_id} = $x;
  }
  return $self->{exon_id};
}


=head3 to_string

Pastes all values of the object to a tab-delimited string.

=cut
sub to_string{
  my $self=shift;
  my $the_string=$self->exon_id."\t".$self->start."\t".$self->end."\t".$self->chromosome_name."\t".$self->strand."\t".$self->coord_system;
  return( $the_string );
}

=head3 print

Prints the values for the Exon.
If an argument "header" is provided it will also print a header.

=cut
sub print{
  my $self=shift;
  my %args = @_;
  if( exists( $args{ header } ) ){
    print "exon_id\tstart\tend\tchromosome_name\tcoord_system\n";
  }
  ## just printing the data... tab separated...
  print $self->to_string."\n";
}


