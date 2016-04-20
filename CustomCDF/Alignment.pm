package CustomCDF::Alignment;
use strict;
use warnings;
use Carp;
use Scalar::Util qw(looks_like_number);
our $VERSION = "0.3.0";

### package variables:



=head1 NAME

Alignment - An alignment of a sequence.

=head1 SYNOPSIS

use CustomCDF::Alignment;

=head1 DESCRIPTION

This represents the alignment of a sequence to a cDNA or chromosome.

Attributes of this object:

=item alignment_id  : the id of the alignment (probes_pk in probe_alignments table).
=item start         : start coordinate(s) of the alignment.
=item end           : the end position(s) of the alignment.
=item probe_id      : the probe_id of the probe.
=item strand        : the strand of the sequence to which the probe is aligned.
=item seq_name      : the name of the sequence, either transcript id or chromosome name.
=item seq_type      : the type of sequence, either cdna or chromosome.
=item mismatches    : number of mismatches of the alignment.


=head2 Methods

=head3 new


=cut

sub new {
  my($class, %args) = @_;
  my $self = bless({}, $class);
  ## my attributes...
  my $start=[];
  if( exists $args{start} ){
    $start = $args{start};
    $self->start($start);
    #croak "Illegal value for start" unless eval { @$start;1 };
    #croak "Illegal value for start" unless ref $start eq "ARRAY";
  }
  my $end=[];
  if( exists $args{end} ){
    $end = $args{end};
    $self->end($end);
    #croak "Illegal value for start" unless eval { @$start;1 };
    #croak "Illegal value for end" unless ref $end eq "ARRAY";
  }
  my $alignment_id = exists $args{alignment_id} ? $args{alignment_id} : "NA";
  my $probe_id = exists $args{probe_id} ? $args{probe_id} : "NA";
  my $strand = exists $args{strand} ? $args{strand} : "NA";
  my $seq_name = exists $args{seq_name} ? $args{seq_name} : "NA";
  my $seq_type = exists $args{seq_type} ? $args{seq_type} : "NA";
  my $mismatches = exists $args{mismatches} ? $args{mismatches} : -1;
  #$self->{start} = $start;
  #$self->{end} = $end;
  $self->{strand} = $strand;
  $self->{alignment_id} = $alignment_id;
  $self->{probe_id} = $probe_id;
  $self->{seq_name} = $seq_name;
  $self->{seq_type} = $seq_type;
  $self->{mismatches} = $mismatches;
  ## get also the reference to the adaptor...if I got fetched from the DB...
  if( exists( $args{adaptor} ) ){
    my $adaptor = $args{adaptor};
    $self->adaptor( $adaptor );
  }
  return $self;
}

=head3 start

my $start = $alignment->start;
$alignment->start( [ 1 ]);
Gets and sets the alignment start position(s). The function returns ALWAYS an array, even if the alignment is un-gapped, i.e. has only a single start position. Also, for setting the start position(s) it is required to submit a reference to an array (either with \@array_name or using [ 1, 2 ,3 ]).

=cut

sub start{
  my $self = shift;
  if( @_ ) {
    my $start = shift;
    croak "Error: illegal value for start" unless ref $start eq 'ARRAY';
    $self->{start} = $start;
  }
  return @{$self->{start}||[]};
}

=head3 end

my $end = $alignment->end;
$alignment->end( [ 1, 2 ] );
Gets and sets the alignment end. The function returns ALWAYS an array, even if the alignment
is un-gapped, i.e. has only a single end position. Also, for setting the end position(s)
it is required to submit a reference to an array (either with \@array_name or using [ 1, 2 ,3 ]).

=cut

sub end{
  my $self = shift;
  if( @_ ) {
    my $end = shift;
    croak "Error: illegal value for start" unless ref $end eq 'ARRAY';
    $self->{end} = $end;
  }
  return @{$self->{end}||[]};
}

=head3 alignment_id

Gets and sets the alignment id.

=cut

sub alignment_id{
  my $self = shift;
  if( @_ ) {
    my $id = shift;
    $self->{alignment_id} = $id;
  }
  return $self->{alignment_id};
}

=head3 probe_id

Gets and sets the probe id.

=cut

sub probe_id{
  my $self = shift;
  if( @_ ) {
    my $id = shift;
    $self->{probe_id} = $id;
  }
  return $self->{probe_id};
}


=head3 strand

Gets and sets the strand. Note: this will return (numeric) +1 and -1 for forward and reverse strand.

=cut

sub strand{
  my $self = shift;
  if( @_ ) {
    my $strand = shift;
    $self->{strand} = $strand;
  }
  ## make sure we always return a number! either +1 or -1!
  my $self_strand=$self->{strand};
  if( looks_like_number( $self_strand ) ){
    if( $self_strand < 0 ){
      return -1;
    }else{
      return 1;
    }
  }else{
    if( $self_strand eq "+" | $self_strand eq "+1" | $self_strand eq "1" | $self_strand eq "fw" | $self_strand eq "forward" ){
      return 1;
    }else{
      return -1;
    }
  }
}

=head3 seq_name

Gets and sets the name of the sequence to which the alignment has been performed.

=cut

sub seq_name{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{seq_name} = $x;
  }
  return $self->{seq_name};
}

=head3 seq_type

Gets and sets the type of the sequence to which the alignment has been performed (chromosome or cdna).

=cut

sub seq_type{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{seq_type} = $x;
  }
  return $self->{seq_type};
}

=head3 mismatches

Gets and sets the number of alignment mismatches.

=cut

sub mismatches{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{mismatches} = $x;
  }
  return $self->{mismatches};
}


#################################################
## more sophisticated calls...


=head3 is_gapped

$alignment->is_gapped();

Returns 1 if the alignment has more than one start position, i.e. is a gapped alignment.

=cut

sub is_gapped{
  my $self = shift;
  if( scalar( $self->start ) > 1 ){
    return 1;
  }else{
    return 0;
  }
}

=head3 adaptor

my $adaptor = $probe->adaptor;
Gets and sets the CustomCDF::ProbeAdaptor adaptor.

=cut
sub adaptor{
  my $self = shift;
  ## this ensures that we do call that only on an object.
  unless (ref $self){
    croak( "Can call adaptor only on an object, not a class!" );
  }
  if( @_ ) {
    my $adaptor = shift;
    if( !ref( $adaptor ) || !$adaptor->isa( 'CustomCDF::ProbeAdaptor' ) ){
      croak( "adaptor has to be a CustomCDF::ProbeAdaptor!" );
    }
    $self->{adaptor} = $adaptor;
  }
  return $self->{adaptor};
}

=head3 have_adaptor

Checks if the object has an adaptor.

=cut
sub have_adaptor{
  my $self = shift;
  return( exists( $self->{adaptor} ) );
}


=head3 to_string

Concatenates the values into a tab separated string.

=cut
sub to_string{
  my $self = shift;
  unless (ref $self){
    croak( "Can call adaptor only on an object, not a class!" )
  }
  my $the_string = $self->alignment_id."\t".$self->probe_id."\t".$self->seq_name."\t".$self->seq_type."\t".join( ",", $self->start)."\t".join( ",", $self->end )."\t".$self->strand."\t".$self->mismatches."\t".$self->is_gapped();
  return( $the_string );
}


=head3 print

Prints the values for the alignment.
If an argument "header" is provided it will also print a header.

=cut
sub print{
  my $self = shift;
  my %args = @_;
  if( exists( $args{ header } ) ){
    print "alignment_id\tprobe_id\tseq_name\tseq_type\tstart\tend\tseq_strand\tmismatches\tgapped\n";
  }
  ## just printing the data... tab separated...
  print $self->to_string."\n";
}


sub DESTROY{
  ## that's called each time we destroy an instance of Probe.
#  print "byebye\n";
}


=head1 AUTHOR

Johannes Rainer <johannes.rainer@i-med.ac.at>

=cut

1;

