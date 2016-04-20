package CustomCDF::Feature;
use strict;
use warnings;
use Carp;
use Scalar::Util qw(looks_like_number);
our $VERSION = "0.3.0";

sub new{
  my($class, %args) = @_;
  my $self = bless({}, $class);
  $self->_init( %args );
  return $self;
}

sub _init{
  my ($self, %args)=@_;
  my $start = exists $args{start} ? $args{start} : -1;
  my $end = exists $args{end} ? $args{end} : -1;
  my $strand = exists $args{strand} ? $args{strand} : 0;
  my $chromosome_name = exists $args{chromosome_name} ? $args{ chromosome_name } : "NA";
  my $coord_system = exists $args{coord_system} ? $args{ coord_system } : "NA";
  if( exists $args{ adaptor } ){
    $self->adaptor( $args{ adaptor } );
  }
  $self->start( $start );
  $self->end( $end );
  $self->strand( $strand );
  $self->chromosome_name( $chromosome_name );
  $self->coord_system( $coord_system );
}

=head3 start

Getter/setter for the (chromosomal) start position of the feature.

=cut
sub start{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{start} = $x;
  }
  return $self->{start};
}


=head3 end

Getter/setter for the (chromosomal) end position of the feature.

=cut
sub end{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{end} = $x;
  }
  return $self->{end};
}


=head3 chromosome_name

Getter/setter for the name of the sequence on which the gene is encoded; in case of a chromosome, the chromosome name (in which case the attribute "coord_system" will be chromosome).

=cut
sub chromosome_name{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{chromosome_name} = $x;
  }
  return $self->{chromosome_name};
}

=head3 strand

Getter/setter for the (chromosomal) strand on which the feature is encoded.

=cut
sub strand{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    my $strand=0;
    ## make sure we save whatever we get as -1, +1:
    if( looks_like_number( $x ) ){
      if( $x < 0 ){
	$strand=-1;
      }else{
	$strand=+1;
      }
    }else{
      if( $x eq "+" | $x eq "+1" | $x eq "1" | $x eq "fw" | $x eq "forward" ){
	$strand=+1;
      }else{
	$strand=-1;
      }
    }
    $self->{strand} = $strand;
  }
  return $self->{strand};
}

=head3 coord_system

Getter/setter for the coordinate system in which the start/end positions are given. In most cases this will be "chromosome", but can also be something else (in which case the chromosome_name will be the name of the contig).

=cut
sub coord_system{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{coord_system} = $x;
  }
  return $self->{coord_system};
}

=head3 adaptor

my $adaptor = $probe->adaptor;
Gets and sets the CustomCDF::GeneAdaptor adaptor.

=cut
sub adaptor{
  my $self = shift;
  ## this ensures that we do call that only on an object.
  unless (ref $self){
    croak( "Can call adaptor only on an object, not a class!" )
  }
  if( @_ ) {
    my $adaptor = shift;
    if( !ref( $adaptor ) || !$adaptor->isa( 'CustomCDF::GeneAdaptor' ) ){
      croak( "adaptor has to be a CustomCDF::GeneAdaptor!" );
    }
    $self->{adaptor} = $adaptor;
  }
  return $self->{adaptor};
}


=head3 length

Gets the length of the feature, i.e. the difference between the end and the start coordinate.

=cut
sub length{
  my $self = shift;
  ## this ensures that we do call that only on an object.
  unless (ref $self){
    croak( "Can call length only on an object, not a class!" )
  }
  my $length = $self->end-$self->start+1;
  return $length;
}

