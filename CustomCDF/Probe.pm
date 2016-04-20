package CustomCDF::Probe;
use strict;
use warnings;
use CustomCDF::Feature;
use Carp;
our $VERSION = "0.3.0";

#our @ISA = qw( CustomCDF::Feature );

### package variables:



=head1 NAME

Probe - A representation of a probe respective its alignment to the genome/cDNA

=head1 SYNOPSIS

use CustomCDF::Probe;
my $probe = CustomCDF::Probe->new();
print $probe->id;

Load a probe from an alignment and annotation database.

my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost",
						  username=>"anonuser",
						  password=>"",
						  dbname=>"homo_sapiens_hugene_cdna_74"
				     );
my $probe_adaptor = $db_adaptor->get_probe_adaptor();
my $probe=$probe_adaptor->fetch_probe( probe_id=>"156902" );
$probe->print();

=head1 DESCRIPTION

This represents a probe from a microarray.

Attributes of this object:

=item id                 : the probe ID.

=item sequence           : the probe's sequence.

=item index_x            : the x-coordinate of the probe on the microarray.

=item index_y            : the y-coordinate of the probe on the microarray.

=item platform           : the microarray type (e.g. HuGene_1.0st).

=item nr_chrom_map       : the number of complete genomic alignments.

=item nr_chrom_map_mmall : the number of genomic alignments with one mismatch.

=item nr_gene_map        : the number of alignments to different genes.

=item nr_gene_map_mm     : the number of alignments to genes with a mismatch.

=item gc_count           : the number of G-C nucleotides in the probe's sequence.

=item cel_index          : the position of the probe within a cel file.

=item adaptor            : the CustomCDF::ProbeAdaptor used to fetch the object from the database.

=head2 Methods

=head3 new

my $probe = Probe->new();
my $hello = Probe->new( id => "probe1" );
Instantiates a Probe object.

=cut

sub new {
  my($class, %args) = @_;
  my $self = bless({}, $class);
  ## my attributes...
  my $id = exists $args{id} ? $args{id} : "NA";
  my $sequence = exists $args{sequence} ? $args{sequence} : "";
  my $x = exists $args{index_x} ? $args{index_x} : 0;
  my $y = exists $args{index_y} ? $args{index_y} : 0;
  my $platform = exists $args{platform} ? $args{platform} : "NA";
  my $nr_chrom_map = exists $args{nr_chrom_map} ? $args{nr_chrom_map} : 0;
  my $nr_chrom_map_mmall = exists $args{nr_chrom_map_mmall} ? $args{nr_chrom_map_mmall} : 0;
  my $nr_gene_map = exists $args{nr_gene_map} ? $args{nr_gene_map} : 0;
  my $nr_gene_map_mm = exists $args{nr_gene_map_mm} ? $args{nr_gene_map_mm} : 0;
  my $gc_count = exists $args{gc_count} ? $args{gc_count} : 0;
  my $cel_index = exists $args{cel_index} ? $args{cel_index} : 0;
  my $alignments=[];
  if( exists $args{alignments} ){
    $alignments=$args{alignments};
    $self->alignments( $alignments );
  }
  ## probe can also contain alignments. that should be a reference to Alignment.
  $self->{id} = $id;
  $self->{sequence} = $sequence;
  $self->{index_x} = $x;
  $self->{index_y} = $y;
  $self->{platform} = $platform;
  $self->{nr_chrom_map} = $nr_chrom_map;
  $self->{nr_chrom_map_mmall} = $nr_chrom_map_mmall;
  $self->{nr_gene_map} = $nr_gene_map;
  $self->{nr_gene_map_mm} = $nr_gene_map_mm;
  $self->{gc_count} = $gc_count;
  $self->{cel_index} = $cel_index;
  ## get also the reference to the adaptor...if I got fetched from the DB...
  if( exists( $args{adaptor} ) ){
    my $adaptor = $args{adaptor};
    $self->adaptor( $adaptor );
  }
  return $self;
}

=head3 id

my $id = $probe->id;
$probe->id("my new id");
Gets and sets the probe id.

=cut

sub id{
  my $self = shift;
  ## this ensures that we do call that only on an object.
  unless (ref $self){
    croak( "Can call id only on an object, not a class!" )
  }
  if( @_ ) {
    my $id = shift;
    $self->{id} = $id;
  }
  return $self->{id};
}

=head3 sequence

my $sequence = $probe->sequence;
$probe->sequence("my new sequence");
Gets and sets the probe sequence.

=cut

sub sequence{
  my $self = shift;
  if( @_ ) {
    my $sequence = shift;
    $self->{sequence} = $sequence;
  }
  return $self->{sequence};
}

=head3 index_x

my $x = $probe->x;
$probe->x(13);
Gets and sets the probe's x-coordinate on the microarray.

=cut

sub index_x{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{index_x} = $x;
  }
  return $self->{index_x};
}

=head3 index_y

my $y = $probe->y;
$probe->y(13);
Gets and sets the probe's y-coordinate on the microarray.

=cut

sub index_y{
  my $self = shift;
  if( @_ ) {
    my $y = shift;
    $self->{index_y} = $y;
  }
  return $self->{index_y};
}

=head3 platform

my $platform = $probe->platform;
$probe->platform("HG-U133_Plus2.0");
Gets and sets the microarray platform.

=cut
sub platform{
  my $self = shift;
  if( @_ ) {
    my $platform = shift;
    $self->{platform} = $platform;
  }
  return $self->{platform};
}


=head3 nr_chrom_map

Getter/setter for the number of (perfect) chromosomal alignments of the probe.

=cut
sub nr_chrom_map{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{nr_chrom_map} = $x;
  }
  return $self->{nr_chrom_map};
}


=head3 nr_chrom_map_mmall

Getter/setter for the number of chromosomal alignments of the probe with one ore more mismatches.

=cut
sub nr_chrom_map_mmall{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{nr_chrom_map_mmall} = $x;
  }
  return $self->{nr_chrom_map_mmall};
}


=head3 nr_gene_map

Getter/setter for the number of genes to which the probe can be aligned.

=cut
sub nr_gene_map{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{nr_gene_map} = $x;
  }
  return $self->{nr_gene_map};
}

=head3 nr_gene_map_mm

Getter/setter for the number of genes to which the probe can be aligned with one (or more) mismatch(es).

=cut
sub nr_gene_map_mm{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{nr_gene_map_mm} = $x;
  }
  return $self->{nr_gene_map_mm};
}


=head3 gc_count

Getter/setter for the G-C content of the probe's sequence.

=cut
sub gc_count{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{gc_count} = $x;
  }
  return $self->{gc_count};
}

=head3 cel_index

Getter/setter for the position of the probe's signal within CEL files.

=cut
sub cel_index{
  my $self = shift;
  if( @_ ) {
    my $x = shift;
    $self->{cel_index} = $x;
  }
  return $self->{cel_index};
}


=head3 alignments

my $algn = CustomCDF::Alignment->new( start=>[1],
				      end=>[10],
				      seq_name=>"test"
				    );
my @alignments = ( $algn );
$test2->alignments( \@alignments );

Gets and sets the alignment(s) for the probe. Alignments has to be an array of CustomCDF::Alignment object. This function will throw and exception if a wrong type of classes is submitted.
Note: this returns alignments that have already been loaded into the object (e.g. with the load_alignments method). In contrast, the get_alignments method will fetch alignments directly from the database (based on specified filters) and will return them directly; without storing them into the object.

=cut

sub alignments{
  my $self = shift;
  if( @_ ) {
    my $algs = shift;
    croak "Error: illegal value for alignments; should be a reference to an array!" unless ref $algs eq "ARRAY";
    my @algs_array = @{ $algs };
    foreach my $the_alg ( @algs_array ){
      #print " ".$the_alg."\n";
      if( !eval { $the_alg->isa( 'CustomCDF::Alignment' ) } ){
	croak "Error: alignments should be a reference to CustomCDF::Alignment objects!";
      }
    }
    $self->{alignments}=$algs;
  }
  return @{$self->{alignments}||[]};
}

=head3 adaptor

my $adaptor = $probe->adaptor;
Gets and sets the CustomCDF::ProbeAdaptor adaptor.

=cut
sub adaptor{
  my $self = shift;
  ## this ensures that we do call that only on an object.
  unless (ref $self){
    croak( "Can call adaptor only on an object, not a class!" )
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


#################################################
## more sophisticated calls...

=head3 load_alignments

This will fetch all alignments (based on some additional filters) of that probe from the database and will store them into the alignments attribute of the class.
Note:this will only query the "probe_alignments" database table unless the option "genomic" is provided.

Optional filters are possible for:

=item seq_name   : e.g. seq_name=>"='ENST00000571914'"

=item seq_type   : e.g. seq_type=>"='cdna'"

=item start      : e.g. start=>"<10000"

=item end        : e.g. end=>">=4"

=item seq_strand : e.g. seq_strand=>"='+'"

=item mismatches : e.g. mismatches=>"=0": defaults to mismatches=>"=0".

=item genomic    (optional)  : fetch all (complete) alignments in genomics coordinates; calls
the fetch_genomic_alignments_for_probe method. This option takes priority over all other options above.

=cut
sub load_alignments{
  ## will call the fetch_alignments_for_probe in the ProbeAdaptor module...
  my $self = shift;
  my %args=@_;
  if( !$self->have_adaptor ){
    croak( "Error: can not load Alignments from the database without a database adaptor!" );
  }
  my @alignments = $self->fetch_alignments( %args );
  $self->alignments( \@alignments );
  return $self;
  # my $adaptor=$self->adaptor;
  # my $probe_id=$self->id;
  # $args{ probe_id } = $probe_id;
  # $args{ chip_type } = $self->chip_type();
  # my @alignments = $adaptor->fetch_alignments_for_probe_id( %args );
  # $self->alignments( \@alignments );
}


=head3 fetch_alignments

This will fetch and return all alignments (based on some additional filters) of that probe from the database.
Note:this will only query the "probe_alignments" database table unless the option "genomic" is provided.

Optional filters are possible for:

=item seq_name   (optional)  : e.g. seq_name=>"='ENST00000571914'"

=item seq_type   (optional)  : e.g. seq_type=>"='cdna'"

=item start      (optional)  : e.g. start=>"<10000"

=item end        (optional)  : e.g. end=>">=4"

=item seq_strand (optional)  : e.g. seq_strand=>"='+'"

=item mismatches (optional)  : e.g. mismatches=>"=0": defaults to mismatches=>"=0";

=item genomic    (optional)  : fetch all (complete) alignments in genomics coordinates; calls
the fetch_genomic_alignments_for_probe method. This option takes priority over all other options above.

=cut
sub fetch_alignments{
  ## will call the fetch_alignments_for_probe in the ProbeAdaptor module...
  my $self = shift;
  my %args=@_;
  if( !$self->have_adaptor ){
    croak( "Error: can not load Alignments from the database without a database adaptor!" );
  }
  my $adaptor=$self->adaptor;
  my $probe_id=$self->id;
  $args{ probe_id } = $probe_id;
  $args{ chip_type } = $self->chip_type();
  my @alignments = $adaptor->fetch_alignments_for_probe_id( %args );
  return( @alignments );
}


=head3 chip_type

Returns "st" for sense target and "at" for antisense target microarrays. The chip type is guessed from the platform name.

=cut
sub chip_type{
  my $self = shift;
  my $platform = $self->platform();
  my $platformlower = "\L$platform";
  my $chiptype="st";
  if( index( $platformlower, "st" ) == -1 ){
    $chiptype="at";
  }
  return $chiptype;
}

=head3 to_string

Concatenates the object's values to a string.

=cut
sub to_string{
  my $self=shift;
  my $the_string=$self->id."\t".$self->index_x."\t".$self->index_y."\t".$self->platform."\t".$self->sequence."\t".$self->nr_chrom_map."\t".$self->nr_chrom_map_mmall."\t".$self->nr_gene_map."\t".$self->nr_gene_map_mm."\t".$self->gc_count."\t".$self->cel_index."\t".scalar( $self->alignments );
  return( $the_string );
}

=head3 print

Prints the values for the probe.
Arguments:

=item header    : print a header.

=item alignment : print also the alignment(s).

=cut
sub print{
  my $self=shift;
  my %args = @_;
  if( exists( $args{ header } ) ){
    print "probe_id\tx\ty\tplatform\tsequence\tnr_chrom_map\tnr_chrom_map_mmall\tnr_gene_map\tnr_gene_map_mm\tgc_count\tcel_index\tloaded_alignments\n";
  }
  ## just printing the data... tab separated...
  print $self->to_string."\n";
  if( exists( $args{ alignment } ) ){
    ## print the alignments if we have some...
    my @algs = $self->alignments();
    if( scalar( @algs ) > 0 ){
      my $do_header=1;
      foreach my $alg ( @algs ){
	if( $do_header==1 ){
	  $do_header=0;
	  $alg->print( header=>1 );
	}else{
	  $alg->print();
	}
      }
    }
  }
}

sub DESTROY{
  ## that's called each time we destroy an instance of Probe.
#  print "byebye\n";
}


=head1 AUTHOR

Johannes Rainer <johannes.rainer@i-med.ac.at>

=cut

1;

