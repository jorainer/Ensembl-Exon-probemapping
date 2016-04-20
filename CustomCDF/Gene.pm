package CustomCDF::Gene;
use strict;
use warnings;
use Carp;
use CustomCDF::Feature;
our $VERSION = "0.3.0";

our @ISA = qw( CustomCDF::Feature );

=head1 NAME

Gene - A gene representation

=head1 SYNOPSIS

This is a representation of a gene. A Gene can be retrieved from an alignment and annotation database using:

my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost",
						  username=>"anonuser",
						  password=>"",
						  dbname=>"homo_sapiens_hugene_cdna_74"
				     );
my $gene_adaptor = $db_adaptor->get_gene_adaptor();
my $gene = $gene_adaptor->fetch_gene( gene_id=>"ENSG00000000003" );

=head1 DESCRIPTION

Is a representation of a gene.
Attributes are:

=item gene_id      : the (Ensembl) gene id

=item gene_name    : the gene name

=item gene_biotype : the biotype

=item coord_system : the coordinate system of the gene; for genes encoded on chromosomes: chromosome

=item transcripts  : an array with transcripts

=item adaptor      : the CustomCDF::GeneAdaptor

CustomCDF::Gene extends CustomCDF::Feature

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
  my $gene_id = exists $args{gene_id} ?  $args{gene_id} : "NA";
  my $gene_name = exists $args{gene_name} ?  $args{gene_name} : "NA";
  my $gene_biotype = exists $args{gene_biotype} ?  $args{gene_biotype} : "NA";
  my $transcripts = [];
  if( exists $args{transcripts} ){
    $transcripts=$args{transcripts};
    $self->transcripts( $transcripts );
  }
  $self->gene_id( $gene_id );
  $self->gene_name( $gene_name );
  $self->gene_biotype( $gene_biotype );
  ## things to be passed to Feature:
  ## coord_system, chromosome_name, strand, adaptor.
  my $strand = exists $args{ transcript_chrom_strand } ? $args{ transcript_chrom_strand } : 0;
  if( !exists $args{ strand } ){
    $args{ strand } = $strand;
  }
  $self->SUPER::_init( %args );
}


=head3 gene_id

Getter/setter for gene id.

=cut
sub gene_id{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{gene_id} = $x;
  }
  return $self->{gene_id};
}

=head3 gene_name

Getter/setter for the gene name.

=cut
sub gene_name{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{gene_name} = $x;
  }
  return $self->{gene_name};
}

=head3 gene_biotype

Getter/setter for the gene biotype.

=cut
sub gene_biotype{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{gene_biotype} = $x;
  }
  return $self->{gene_biotype};
}


=head3 transcripts

Getter/setter for the transcripts attribute of the gene. Note that transcripts can be either loaded directly when the gene is fetched from the database, or later using the load_transcripts method. The fetch_transcripts method on the other hand will also retrieve the transcripts for the gene, will however not store them into the transcripts attribute.

=cut
sub transcripts{
  my $self = shift;
  if( @_ ) {
    my $transcripts = shift;
    croak "Error: illegal value for transcripts; should be a reference to an array!" unless ref $transcripts eq "ARRAY";
    my @transcript_array = @{ $transcripts };
    foreach my $the_transcript ( @transcript_array ){
      if( !eval { $the_transcript->isa( 'CustomCDF::Transcript' ) } ){
	croak "Error: transcripts should be a reference to CustomCDF::Transcript objects!";
      }
    }
    $self->{transcripts}=$transcripts;
  }
  return @{$self->{transcripts}||[]};
}


=head3 fetch_transcripts

Fetch all transcripts for that gene from the database. Each call of this method will result in a new database query.
Parameters:

=item load_exons   (optional)  : optionally load also all exons for each transcript.

=cut
sub fetch_transcripts{
  my $self = shift;
  my %args = @_;
  my $load_exons = exists( $args{ load_exons } );
  if( $load_exons ){
    return ( $self->adaptor()->fetch_transcripts_for_gene( gene_id=>$self->gene_id, load_exons=>1 ) );
  }else{
    return ( $self->adaptor()->fetch_transcripts_for_gene( gene_id=>$self->gene_id ) );
  }
}

=head3 load_transcripts

Fetch all transcripts for that gene from the database and store them to the transcripts attribute. Each call of this method will result in a new database query.
Parameters:

=item load_exons   (optional)  : optionally load also all exons for each transcript.

=cut
sub load_transcripts{
  my $self = shift;
  my %args = @_;
  my @transcripts = $self->fetch_transcripts(%args);
  $self->transcripts( \@transcripts );
  return @transcripts;
}



=head3 fetch_exons

Retrieves all (unique) exons for the gene. Returns an array of CustomCDF::Exon, exons are ordered by their chromosomal start position.

=cut
sub fetch_exons{
  my $self = shift;
  return ( $self->adaptor()->fetch_exons_for_gene( gene_id=>$self->gene_id ) );
}

=head3 to_string

Pastes all values of the object to a tab-delimited string.

=cut
sub to_string{
  my $self=shift;
  my $have_transcripts = scalar $self->transcripts();
  my $the_string=$self->gene_id."\t".$self->gene_name."\t".$self->gene_biotype."\t".$self->chromosome_name."\t".$self->strand."\t".$self->coord_system;
  if( $have_transcripts > 0 ){
    $the_string=$the_string."\t".$have_transcripts;
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
    my $header = "gene_id\tgene_name\tgene_biotype\tchromosome_name\tstrand\tcoord_system";
    my $have_transcripts = scalar $self->transcripts();
    if( $have_transcripts > 0 ){
      $header = $header."\tn_transcripts";
    }
    print $header."\n";
  }
  ## just printing the data... tab separated...
  print $self->to_string."\n";
}

1;

