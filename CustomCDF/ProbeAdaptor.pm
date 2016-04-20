package CustomCDF::ProbeAdaptor;
use CustomCDF::Probe;
use CustomCDF::Alignment;
use strict;
use DBI;
use IO::File;
use warnings;
use Carp;
use Config::Simple;
our $VERSION = "0.3.0";

our @ISA = qw( CustomCDF::DBAdaptor );
### package variables:
## the queriies:
my $query_fetch_probe = "select probe_id,sequence,x,y,platform from probe_information where probe_id=?";

## load all probes for exon. NOTE: only probes with complete alignments within the exon are returned! no splice junction probes!
my $query_fetch_all_probes_for_exon = "select probe_information.probe_id, sequence, x, y, platform,nr_chrom_map, nr_chrom_map_mmall, nr_gene_map, nr_gene_map_mm, gc_count, cel_index from (select distinct  probe_id from probe_to_exon_transcript where exon_id=? and is_junction_alignment=0 order by probe_alignment_exon_start) as ptet join probe_information on (ptet.probe_id=probe_information.probe_id);";

## load probes for exon: only probes with one (or less) genomic alignments and no additional alignments with up to one mismatch: NOTE: only probes with complete alignments within the exon are returned! no splice junction probes!
my $query_fetch_probes_for_exon = "select probe_information.probe_id, sequence, x, y, platform,nr_chrom_map, nr_chrom_map_mmall, nr_gene_map, nr_gene_map_mm, gc_count, cel_index from (select distinct probe_id from probe_to_exon_transcript where exon_id=? and is_junction_alignment=0 order by probe_alignment_exon_start) as ptet join probe_information on (ptet.probe_id=probe_information.probe_id) where nr_chrom_map<=1 and nr_chrom_map_mmall=0;";

## in addition to the fetch_probe_for_exon we are also loading the alignment of the probe WITHIN THE EXONS! We are fetching this information from the probe_to_exon_transcript table.
my $query_fetch_probes_with_alignments_for_exon = "select probe_information.probe_id, sequence, x, y, platform,nr_chrom_map, nr_chrom_map_mmall, nr_gene_map, nr_gene_map_mm, gc_count, cel_index, probe_alignment_exon_start, probe_alignment_exon_end, probes_fk, chromosome_strand from (select probe_id, probe_alignment_exon_start, probe_alignment_exon_end, chromosome_strand, probes_fk from probe_to_exon_transcript where exon_id=? and is_junction_alignment=0 order by probe_alignment_exon_start) as ptet join probe_information on (ptet.probe_id=probe_information.probe_id) where nr_chrom_map<=1 and nr_chrom_map_mmall=0;";

## load all probes for a transcript: all probes listed in the probe_to_exon_transcript table for a given transcript are returned and the results is joined with the probe_information table.
my $query_fetch_all_probes_for_transcript = "select probe_information.probe_id, sequence, x, y, platform, nr_chrom_map, nr_chrom_map_mmall, nr_gene_map, nr_gene_map_mm, gc_count, cel_index from (select distinct probe_id from probe_to_exon_transcript where transcript_id=? order by probe_alignment_exon_start) as ptet_at join probe_information on (ptet_at.probe_id=probe_information.probe_id);";

## load "good" probes for a transcript: only probes with one (or less) genomic alignments and no additional alignments with up to one mismatch: NOTE: this query might also return exon junction probes.
my $query_fetch_probes_for_transcript = "select probe_information.probe_id, sequence, x, y, platform, nr_chrom_map, nr_chrom_map_mmall, nr_gene_map, nr_gene_map_mm, gc_count, cel_index from (select distinct probe_id from probe_to_exon_transcript where transcript_id=? order by probe_alignment_exon_start) as ptet_t join probe_information on (ptet_t.probe_id=probe_information.probe_id) where nr_chrom_map<=1 and nr_chrom_map_mmall=0;";

## load all probes for gene. NOTE: only probes with complete alignments within the exon are returned! no splice junction probes!
my $query_fetch_all_probes_for_gene = "select probe_information.probe_id, sequence, x, y, platform,nr_chrom_map, nr_chrom_map_mmall, nr_gene_map, nr_gene_map_mm, gc_count, cel_index from (select distinct probe_id from probe_to_exon_transcript where gene_id=? and is_junction_alignment=0 order by probe_alignment_exon_start) as ptet_ag join probe_information on (ptet_ag.probe_id=probe_information.probe_id);";

## load probes for gene: only probes with one (or less) genomic alignments and no additional alignments with up to one mismatch: NOTE: only probes with complete alignments within the exon are returned! no splice junction probes!
my $query_fetch_probes_for_gene = "select probe_information.probe_id, sequence, x, y, platform, nr_chrom_map, nr_chrom_map_mmall, nr_gene_map, nr_gene_map_mm, gc_count, cel_index from (select distinct probe_id from probe_to_exon_transcript where gene_id=? and is_junction_alignment=0 order by probe_alignment_exon_start) as ptet_g join probe_information on (ptet_g.probe_id=probe_information.probe_id) where nr_chrom_map<=1 and nr_chrom_map_mmall=0;";


## get the exon start/end coordinates.
my $query_fetch_exon_startend = "select distinct exon_id, exon_chrom_start, exon_chrom_end, transcript_chrom_strand, chromosome_name from exon_transcript where exon_id=?";
## get the transcript start/end coordinates.
my $query_fetch_transcript_startend = "select distinct transcript_id, transcript_chrom_start, transcript_chrom_end, transcript_chrom_strand, chromosome_name from exon_transcript where transcript_id=?";

## fetch genomic alignment(s) for a probe; only alignments without mismatches are returned.
my $query_fetch_genomic_alignments_for_probe = "select probes_pk, seq_name, seq_type, start, end, seq_strand from probe_alignments where probe_id=? and missmatches=0 and seq_type='chromosome'";

## fetch all we need to calculate genomic alignments for alignments against cdna.
my $query_fetch_for_genomic_alignments_for_probe = "select distinct chromosome_strand, probe_id, probe_alignment_exon_start, probe_alignment_exon_end, exon_chrom_start, exon_chrom_end, chromosome_name, probes_fk from (select probes_fk, chromosome_strand, probe_id, exon_id, probe_alignment_exon_start, probe_alignment_exon_end from probe_to_exon_transcript where probe_id=?) as ptet2 join exon_transcript on ( ptet2.exon_id=exon_transcript.exon_id );";

=head1 NAME

ProbeAdaptor - An object to retrieve probes from the database.

=head1 SYNOPSIS

use CustomCDF::ProbeAdaptor;
use CustomCDF::DBAdaptor;
my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost",
						  username=>"anonuser",
						  password=>"",
						  dbname=>"homo_sapiens_hugene_cdna_74"
				     );
my $probe_adaptor = $db_adaptor->get_probe_adaptor();

my $probe_adaptor = CustomCDF::ProbeAdaptor->new();


=head1 DESCRIPTION

Get Probe objects from the custom CDF annotation database.

=head2 Methods

=head3 new

The preferred way to generate a ProbeAdaptor is via the CustomCDF::DBAdaptor.

Alternatively, a ProbeAdaptor can be generated using the new method.

my $probe = CustomCDF::ProbeAdaptor->new( host=>"localhost", username=>"anonuser", password=>"", dbname=>"homo_sapiens_hugene_cdna_73");

Instantiates a ProbeAdaptor object. Alternatively a configuration file can be submitted using the config_file argument.

=cut

sub new {
  my $class = shift;
  #my($class, %args) = @_;
  my $self = bless({}, $class);
  ## according to http://www.perlmonks.org/?node_id=176963 it's better to separate new from initialization.
  ## call new on DBAdaptor
  #$self->SUPER::new(@_);
  ## initialize all
  $self->_init(@_);
  return $self;
}

sub _init{
  my ($self, %args)=@_;
  $self->SUPER::_init(%args);
  ## and now prepare all cached queries.
  my $dbcon=$self->dbcon();
  my $pq_fetch_probe=$dbcon->prepare_cached( $query_fetch_probe ) or croak $dbcon->errstr;
  $self->{pq_fetch_probe}=$pq_fetch_probe;
  ## exon related stuff...
  my $pq_fetch_all_probes_for_exon=$dbcon->prepare_cached( $query_fetch_all_probes_for_exon ) or croak $dbcon->errstr;
  $self->{pq_fetch_all_probes_for_exon}=$pq_fetch_all_probes_for_exon;
  my $pq_fetch_probes_for_exon=$dbcon->prepare_cached( $query_fetch_probes_for_exon ) or croak $dbcon->errstr;
  $self->{pq_fetch_probes_for_exon}=$pq_fetch_probes_for_exon;
  my $pq_fetch_probes_with_alignments_for_exon=$dbcon->prepare_cached( $query_fetch_probes_with_alignments_for_exon ) or croak $dbcon->errstr;
  $self->{pq_fetch_probes_with_alignments_for_exon}=$pq_fetch_probes_with_alignments_for_exon;
  ## transcript related stuff
  my $pq_fetch_probes_for_transcript=$dbcon->prepare_cached( $query_fetch_probes_for_transcript ) or croak $dbcon->errstr;
  $self->{pq_fetch_probes_for_transcript}=$pq_fetch_probes_for_transcript;
  my $pq_fetch_all_probes_for_transcript=$dbcon->prepare_cached( $query_fetch_all_probes_for_transcript ) or croak $dbcon->errstr;
  $self->{pq_fetch_all_probes_for_transcript}=$pq_fetch_all_probes_for_transcript;
  ## gene related stuff
  my $pq_fetch_probes_for_gene=$dbcon->prepare_cached( $query_fetch_probes_for_gene ) or croak $dbcon->errstr;
  $self->{pq_fetch_probes_for_gene}=$pq_fetch_probes_for_gene;
  my $pq_fetch_all_probes_for_gene=$dbcon->prepare_cached( $query_fetch_all_probes_for_gene ) or croak $dbcon->errstr;
  $self->{pq_fetch_all_probes_for_gene}=$pq_fetch_all_probes_for_gene;
  ## query to transform cdna alignment to genomic alignment.
  my $pq_fetch_for_genomic_alignment_for_probe=$dbcon->prepare_cached( $query_fetch_for_genomic_alignments_for_probe ) or croak $dbcon->errstr;
  $self->{pq_fetch_for_genomic_alignment_for_probe}=$pq_fetch_for_genomic_alignment_for_probe;
  my $pq_fetch_genomic_alignment_for_probe=$dbcon->prepare_cached( $query_fetch_genomic_alignments_for_probe ) or croak $dbcon->errstr;
  $self->{pq_fetch_genomic_alignment_for_probe}=$pq_fetch_genomic_alignment_for_probe;
}

#####
## getters for the cached queries.
sub _pq_fetch_probe{
  my $self = shift;
  return( $self->{pq_fetch_probe} );
}
sub _pq_fetch_all_probes_for_exon{
  my $self = shift;
  return( $self->{pq_fetch_all_probes_for_exon } );
}
sub _pq_fetch_probes_for_exon{
  my $self = shift;
  return( $self->{pq_fetch_probes_for_exon } );
}
sub _pq_fetch_probes_for_transcript{
  my $self = shift;
  return( $self->{pq_fetch_probes_for_transcript } );
}
sub _pq_fetch_all_probes_for_transcript{
  my $self = shift;
  return( $self->{pq_fetch_all_probes_for_transcript } );
}
sub _pq_fetch_probes_for_gene{
  my $self = shift;
  return( $self->{pq_fetch_probes_for_gene } );
}
sub _pq_fetch_all_probes_for_gene{
  my $self = shift;
  return( $self->{pq_fetch_all_probes_for_gene } );
}
sub _pq_fetch_probes_with_alignments_for_exon{
  my $self = shift;
  return( $self->{pq_fetch_probes_with_alignments_for_exon } );
}
sub _pq_fetch_for_genomic_alignment_for_probe{
  my $self = shift;
  return( $self->{pq_fetch_for_genomic_alignment_for_probe } );
}
sub _pq_fetch_genomic_alignment_for_probe{
  my $self = shift;
  return( $self->{pq_fetch_genomic_alignment_for_probe } );
}

sub DESTROY{
  ## that's called each time we destroy an instance of Probe.
#  print "byebye\n";
  my $self = shift;
  ## finish all prepared queries.
  if( defined $self->dbcon ){
    $self->_pq_fetch_probe->finish;
    ## exon related:
    $self->_pq_fetch_all_probes_for_exon->finish;
    $self->_pq_fetch_probes_for_exon->finish;
    $self->_pq_fetch_probes_with_alignments_for_exon->finish;
    $self->_pq_fetch_for_genomic_alignment_for_probe->finish;
    $self->_pq_fetch_genomic_alignment_for_probe->finish;
    $self->_pq_fetch_all_probes_for_transcript->finish;
    $self->_pq_fetch_probes_for_transcript->finish;
    $self->_pq_fetch_all_probes_for_gene->finish;
    $self->_pq_fetch_probes_for_gene->finish;
  }
}



################################################################
##   Probe related stuff.
##
################################################################

=head3 fetch_probe


Retrieves a Probe object from the database for the probe '123'.
Internally we are fetching the probe from the probe_information table.
Parameters:

=item probe_id    (mandatory)  : the id of the probe.

=item load_alignments          (optional)  : fetch and load all alignments for the probe as reported in the "probe_alignments" database table. See options below for further filters.

=item seq_name                 (optional)  : e.g. seq_name=>"='ENST00000571914'"

=item seq_type                 (optional)  : e.g. seq_type=>"='cdna'"

=item start                    (optional)  : e.g. start=>"<10000"

=item end                      (optional)  : e.g. end=>">=4"

=item seq_strand               (optional)  : e.g. seq_strand=>"='+'"

=item mismatches               (optional)  : e.g. mismatches=>"=0"; defaults to mismatches=>"=0";

=item load_genomic_alignments  (optional)  : fetch all (complete) alignments in genomics coordinates; calls the fetch_genomic_alignments_for_probe method. The "load_alignments" option is preferred if both "load_alignments" and "load_genomic_alignments" are called.

=item chip_type  (optional)  : the chip type, either st or at; see help for fetch_genomic_alignments_for_probe.


Example:
my $probe = $probe_adaptor->fetch_probe( probe_id=>'123');

=cut
sub fetch_probe{
  my $self = shift;
  ## this ensures that we do call that only on an object.
  unless (ref $self){
    croak( "Can call fetch_by_probe_id only on an object, not a class!" )
  }
  my %args = @_;
  croak "Error: the probe id has to be speficied with probe_id=>!" unless exists( $args{ probe_id } );
  my $probe_id = $args{probe_id};
  my $load_alignments = exists( $args{ load_alignments } );
  my $load_genomic_alignments = exists( $args{ load_genomic_alignments } );
  my $dbcon = $self->dbcon;
  $self->_pq_fetch_probe->execute($probe_id) or croak $dbcon->errstr;
  my %rowhash = %{ $self->_pq_fetch_probe->fetchrow_hashref };
  if( %rowhash ){
    $rowhash{ id } = $rowhash{ probe_id };
    my $probe = CustomCDF::Probe->new( %rowhash, adaptor=>$self );
    if( $load_alignments ){
      ## load alignments using the fetch_alignments_for_probe_id
      my @alignments = $self->fetch_alignments_for_probe_id( %args );
      $probe->alignments( \@alignments );
    }
    if( $load_genomic_alignments ){
      ## require the chip type and call fetch_genomic_alignments_for_probe
      $args{ chip_type } = $probe->chip_type();
      my @alignments = $self->fetch_genomic_alignments_for_probe( %args );
      $probe->alignments( \@alignments );
    }
    return $probe;
  }
}

=head3 fetch_all_probes

Fetch all probes from the probe_information table optionally matching one or more settings.

Optional filters:

=item nr_gene_map       (optional) : e.g. nr_gene_map=>">=1"

=item nr_chrom_map      (optional) : e.g. nr_chrom_map=>"=1"

=item nr_chom_map_mmall (optional) : e.g. nr_chrom_map_mmall=>"=0"

=item nr_gene_map_mm    (optional) : e.g. nr_gene_map_mm=>"=0"

=item gc_count          (optional) : e.g. gc_count=>"<19"

Examples:

use CustomCDF::ProbeAdaptor;
use CustomCDF::Probe;
## directly create the adaptor.
my $probe_adaptor = CustomCDF::ProbeAdaptor->new( host=>"localhost",
						  username=>"anonuser",
						  password=>"",
						  dbname=>"homo_sapiens_hugene_cdna_74"
				     );

## preferred way:
my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost",
						  username=>"anonuser",
						  password=>"",
						  dbname=>"homo_sapiens_hugene_cdna_74"
				     );
my $probe_adaptor = $db_adaptor->get_probe_adaptor();

## fetching all probes
my @probes = $probe_adaptor->fetch_all_probes();
print "got: ".scalar(@probes)." probes\n";


## fetching probes which map to a single gene and do not have more than 1 alignment to the genome (also considering a missmatch)
my @probes = $probe_adaptor->fetch_all_probes( nr_gene_map=>"=1", nr_gene_map_mm=>"<=0", nr_chrom_map=>"<=1", nr_chrom_map_mmall=>"=0" );
print "got: ".scalar(@probes)." probes\n";

=cut
sub fetch_all_probes{
  my $self = shift;
#  my($class, %args) = @_;
  my %args = @_;
  my $nr_gene_map = exists $args{nr_gene_map} ? $args{nr_gene_map} : "-1";
  my $nr_chrom_map = exists $args{nr_chrom_map} ? $args{nr_chrom_map} : "-1";
  my $nr_chrom_map_mmall = exists $args{nr_chrom_map_mmall} ? $args{nr_chrom_map_mmall} : "-1";
  my $nr_gene_map_mm = exists $args{nr_gene_map_mm} ? $args{nr_gene_map_mm} : "-1";
  my $gc_count = exists $args{gc_count} ? $args{gc_count} : "-1";
  my $query_string = "select probe_id, x, y, platform, sequence, nr_chrom_map, nr_chrom_map_mm1, nr_chrom_map_mmall, nr_gene_map, nr_gene_map_mm, gc_count, cel_index from probe_information";
  my $condition = " where ";
  if( $nr_gene_map ne "-1" ){
    $query_string = $query_string.$condition."nr_gene_map".$nr_gene_map;
    $condition = " and ";
  }
  if( $nr_gene_map_mm ne "-1" ){
    $query_string = $query_string.$condition."nr_gene_map_mm".$nr_gene_map_mm;
    $condition = " and ";
  }
  if( $nr_chrom_map ne "-1" ){
    $query_string = $query_string.$condition."nr_chrom_map".$nr_chrom_map;
    $condition = " and ";
  }
  if( $nr_chrom_map_mmall ne "-1" ){
    $query_string = $query_string.$condition."nr_chrom_map_mmall".$nr_chrom_map_mmall;
    $condition = " and ";
  }
  if( $gc_count ne "-1" ){
    $query_string = $query_string.$condition."gc_count".$gc_count;
    $condition = " and ";
  }
  $query_string=$query_string.";";
  my $db_con=$self->{dbcon};
  my $query = $db_con->prepare( $query_string ) or croak $db_con->errstr;
  $query->execute() or croak $db_con->errstr;
  my @results = ();
  ## looping through the probes and storing them into the array...
  while ( my $hashref = $query->fetchrow_hashref ){
    my %rowhash = %{$hashref};
    $rowhash{id}=$rowhash{probe_id};
    $rowhash{index_x}=$rowhash{x};
    $rowhash{index_y}=$rowhash{y};
    push( @results, CustomCDF::Probe->new( %rowhash, adaptor=>$self ) );
  }
  # while ( my @row = $query->fetchrow_array){
  #   push( @results, $self->row_to_probe( @row ) );
  # }
  return( @results );
}


=head3 fetch_probes_for_exon

Fetch all probes from the database that allow to target the exon. Note that only "good" probes are returned, i.e. probes that perfectly match the exon and that do not have any additional genomic alignment allowing also up to one mismatch (the "all" parameter allows to retrieve all probes, also "bad" ones). Note that only probes with comlpete alignments within the exon are returned, thus excluding exon-junction probes.
Parameters:

=item exon_id                 (mandatory): the id of the exon.

=item all                     (optional)   : retrieves all probes perfectly aligning within the given exon. Note: this retrieves also probes with multiple genomic alignments, with or without mismatches. Also, this parameter is mutually exclusive with the parameter load_alignments!

=item load_alignments         (optional) : whether alignments should also be loaded. The reported alignments are within the exon, thus the start and end positions are relative to the exon (they are not chromosomal coordinates).

=item load_genomic_alignments (optional) : returns all alignments of the probe in genomic coordinates. Note: this will result in an additional SQL query for each probe. load_alignments if preferred if both load_alignments and load_genomic_alignments are called.

=cut
sub fetch_probes_for_exon{
  ## fetch all probes that match the exon using the probe_to_exon_transcript table.
  ## Note: probe_to_exon_transcript contains only probes that match the gene/transcript/exons without mismatches, however, does also contain potentially probes with multiple alignments. To exclude those we have to join the table with probe_information.
  my $self = shift;
  ## here we assume some additional stuff...
  my %args = @_;
  croak "Error: the exon id has to be speficied with exon_id=>!" unless exists( $args{ exon_id } );
  my $exon_id = $args{exon_id};
  my $load_alignments = exists( $args{ load_alignments } );
  my $load_genomic_alignments = exists( $args{ load_genomic_alignments } );
  my $all = exists( $args{ all } );
  my $dbcon = $self->dbcon;
  ## different possibilities:
  ## load "all" probes (without alignments!).
  my @probes=();
  if( $all ){
    $self->_pq_fetch_all_probes_for_exon->execute( $exon_id ) or croak $dbcon->errstr;
    while ( my $hashref = $self->_pq_fetch_all_probes_for_exon->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      ## doing a little "re-formatting"
      $rowhash{ id } = $rowhash{ probe_id };
      $rowhash{ index_x } = $rowhash{ x };
      $rowhash{ index_y } = $rowhash{ y };
      $rowhash{ adaptor } = $self;
      push( @probes, CustomCDF::Probe->new( %rowhash ) );
    }
    return( @probes );
  }
  ## load "good" probes without alignments.
  if( !$load_alignments && !$load_genomic_alignments ){
    $self->_pq_fetch_probes_for_exon->execute( $exon_id ) or croak $dbcon->errstr;
    while ( my $hashref = $self->_pq_fetch_probes_for_exon->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      ## doing a little "re-formatting"
      $rowhash{ id } = $rowhash{ probe_id };
      $rowhash{ index_x } = $rowhash{ x };
      $rowhash{ index_y } = $rowhash{ y };
      $rowhash{ adaptor } = $self;
      push( @probes, CustomCDF::Probe->new( %rowhash ) );
    }
    return( @probes );
  }
  ## load "good" probes with alignment(s). Note that we can afford doing this with one query, since we expect a 1:1 mapping of probe:alignment.
  if( $load_alignments ){
    $self->_pq_fetch_probes_with_alignments_for_exon->execute( $exon_id ) or croak $dbcon->errstr;
    while ( my $hashref = $self->_pq_fetch_probes_with_alignments_for_exon->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      ## doing a little "re-formatting"
      $rowhash{ id } = $rowhash{ probe_id };
      $rowhash{ index_x } = $rowhash{ x };
      $rowhash{ index_y } = $rowhash{ y };
      $rowhash{ adaptor } = $self;
      my $probe = CustomCDF::Probe->new( %rowhash );
      ## for the alignment: got probe_alignment_exon_start, probe_alignment_exon_end, chromosome_strand, probes_fk.
      ## formatting for the alignment:
      ## define the strand: since we are reporting the alignment WITHIN THE EXON, strand is - for ST and + for AT arrays.
      if( $probe->chip_type eq "st" ){
	$rowhash{ strand } = -1;
      }else{
	$rowhash{ strand } = 1;
      }
      $rowhash{ alignment_id } = $rowhash{ probes_fk };
      $rowhash{ seq_name } = $exon_id;  ## coordinates are relative to the exon
      $rowhash{ seq_type } = "cdna";    ##
      $rowhash{ mismatches } = 0;       ## that's the way to go.
      $rowhash{ start } = [ $rowhash{ probe_alignment_exon_start } ];
      $rowhash{ end } = [ $rowhash{ probe_alignment_exon_end } ];
      my $algn = CustomCDF::Alignment->new( %rowhash );
      $probe->alignments( [ $algn ] );
      push( @probes, $probe );
    }
    return( @probes );
  }
  ## load the genomic alignments for the probe
  if( $load_genomic_alignments ){
    $self->_pq_fetch_probes_for_exon->execute( $exon_id ) or croak $dbcon->errstr;
    while ( my $hashref = $self->_pq_fetch_probes_for_exon->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      ## doing a little "re-formatting"
      $rowhash{ id } = $rowhash{ probe_id };
      $rowhash{ index_x } = $rowhash{ x };
      $rowhash{ index_y } = $rowhash{ y };
      $rowhash{ adaptor } = $self;
      my $probe = CustomCDF::Probe->new( %rowhash );
      my @algs = $self->fetch_genomic_alignments_for_probe( probe_id=>$probe->id(), chip_type=>$probe->chip_type );
      $probe->alignments( \@algs );
      push( @probes, $probe );
    }
    return( @probes );
  }
}


=head3 fetch_probes_for_transcript

Fetch all probes from the database that allow to target the transcript. Note that by default only "good" probes are returned, i.e. probes that perfectly match the transcript and that do not have any additional genomic alignment allowing also up to one mismatch. Parameter "all" allows to retrieve all probes even those with multiple partial or complete alignments.
The probe to transcript information recorded in the "probe_to_exon_transcript" table is derived from the alignment of the probes to the cDNA. Thus probes with alignments over splice junctions, which result in two rows in the "probe_to_exon_transcript" table are always correctly assigned to the respective transcript. In other words, querying for probes of a given transcript in the "probe_to_exon_transcript" table will never return splice junction probes that match only partially to the exons of the transcript.

Parameters:

=item transcript_id            (mandatory): the id of the transcript.

=item all                      (optional) : retrieves also probes that could have additional partial or complete genomic alignments.

=item load_alignments          (optional) : whether alignments should also be loaded (only alignments without mismatches are returned).
Additional parameters to filter alignments can be specified (see fetch_alignments_for_probe_id).

=item load_genomic_alignments  (optional) : whether genomic alignments should be loaded (only alignments without mismatches are returned). Note that this causes an additional SQL query for each probe. Also, this parameter is mutually exclusive with "load_alignments", which is preferred if both are set.

=cut
sub fetch_probes_for_transcript{
  ## do it the "slow" way: first get all probes and call load_alignments on those, depending on the settings above.
  my $self = shift;
  my %args = @_;
  croak "Error: the transcript id has to be speficied with transcript_id=>!" unless exists( $args{ transcript_id } );
  my $transcript_id = $args{transcript_id};
  ## define which probes and what we want to get.
  my $all = exists( $args{ all } );     ## load all(!) probes.
  my $load_alignments = exists( $args{ load_alignments } );
  my $load_genomic_alignments = exists( $args{ load_genomic_alignments } );
  my $dbcon = $self->dbcon;
  ## TODO: query probe_to_exon_transcripts for probe_id, transcript_id, probes_fk, two different queries, depending on all.
  ## for each probe, call fetch_alignments... takes longer, but keeps source cleaner.
  my @probes=();
  if( $all ){
    $self->_pq_fetch_all_probes_for_transcript->execute( $transcript_id ) or croak $dbcon->errstr;
    while ( my $hashref = $self->_pq_fetch_all_probes_for_transcript->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      ## doing a little "re-formatting"
      $rowhash{ id } = $rowhash{ probe_id };
      $rowhash{ index_x } = $rowhash{ x };
      $rowhash{ index_y } = $rowhash{ y };
      $rowhash{ adaptor } = $self;
      my $probe = CustomCDF::Probe->new( %rowhash );
      if( $load_alignments ){
	my @algs = $self->fetch_alignments_for_probe_id( probe_id=>$rowhash{ probe_id } );
	$probe->alignments( \@algs );
      }
      if( $load_genomic_alignments ){
	my @algs = $self->fetch_genomic_alignments_for_probe( probe_id=>$rowhash{ probe_id }, chip_type=>$probe->chip_type() );
	$probe->alignments( \@algs );
      }
      push( @probes, $probe );
    }
  }else{
    $self->_pq_fetch_probes_for_transcript->execute( $transcript_id ) or croak $dbcon->errstr;
    while ( my $hashref = $self->_pq_fetch_probes_for_transcript->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      ## doing a little "re-formatting"
      $rowhash{ id } = $rowhash{ probe_id };
      $rowhash{ index_x } = $rowhash{ x };
      $rowhash{ index_y } = $rowhash{ y };
      $rowhash{ adaptor } = $self;
      my $probe = CustomCDF::Probe->new( %rowhash );
      if( $load_alignments ){
	my %dummyhash = %args;
	$dummyhash{ probe_id } = $rowhash{ probe_id };
	#my @algs = $self->fetch_alignments_for_probe_id( probe_id=>$rowhash{ probe_id } );
	my @algs = $self->fetch_alignments_for_probe_id( %dummyhash );
	$probe->alignments( \@algs );
      }
      if( $load_genomic_alignments ){
	my @algs = $self->fetch_genomic_alignments_for_probe( probe_id=>$rowhash{ probe_id }, chip_type=>$probe->chip_type() );
	$probe->alignments( \@algs );
      }
      push( @probes, $probe );
    }
  }
  return( @probes );
}


=head3 fetch_probes_for_gene

Fetch all probes from the database that allow to target the gene. Note that by default only "good" probes are returned, i.e. probes that perfectly match the gene and that do not have any additional genomic alignment allowing also up to one mismatch. Parameter "all" allows to retrieve all probes even those with multiple partial or complete alignments. Note: this method does not return splice junction probes! These probes are only retrieved for transcripts, not for genes or exons.

Parameters:

=item gene_id                  (mandatory): the id of the gene.

=item all                      (optional) : retrieves also probes that could have additional partial or complete genomic alignments.

=item load_alignments          (optional) : whether alignments should also be loaded (only alignments without mismatches are returned). Note: this causes an additional SQL query for each probe. Also, alignments to cdna and chromosome are returned. Thus alignments of the same probe can be returned twice, once in chromosomal coordinates, once in cDNA coordinates.
Additional parameters to filter alignments can be specified (see fetch_alignments_for_probe_id).

=item load_genomic_alignments  (optional) : whether genomic alignments should be loaded (only alignments without mismatches are returned). Note that this causes an additional SQL query for each probe. Also, this parameter is mutually exclusive with "load_alignments", which is preferred if both are set. For each probe only a single alignment in genomic coordinates will be returned. Genomic coordinates are retrieved directly from the database and are in addition calculated based on alignment to cDNA and the genomic coordinates of the cDNA.

=cut
sub fetch_probes_for_gene{
  ## do it the "slow" way: first get all probes and call load_alignments on those, depending on the settings above.
  my $self = shift;
  my %args = @_;
  croak "Error: the gene id has to be speficied with gene_id=>!" unless exists( $args{ gene_id } );
  my $gene_id = $args{gene_id};
  ## define which probes and what we want to get.
  my $all = exists( $args{ all } );     ## load all(!) probes.
  my $load_alignments = exists( $args{ load_alignments } );
  my $load_genomic_alignments = exists( $args{ load_genomic_alignments } );
  my $dbcon = $self->dbcon;
  ## for each probe, call fetch_alignments... takes longer, but keeps source cleaner.
  my @probes=();
  if( $all ){
    $self->_pq_fetch_all_probes_for_gene->execute( $gene_id ) or croak $dbcon->errstr;
    while ( my $hashref = $self->_pq_fetch_all_probes_for_gene->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      ## doing a little "re-formatting"
      $rowhash{ id } = $rowhash{ probe_id };
      $rowhash{ index_x } = $rowhash{ x };
      $rowhash{ index_y } = $rowhash{ y };
      $rowhash{ adaptor } = $self;
      my $probe = CustomCDF::Probe->new( %rowhash );
      if( $load_alignments ){
	my @algs = $self->fetch_alignments_for_probe_id( probe_id=>$rowhash{ probe_id } );
	$probe->alignments( \@algs );
      }
      if( $load_genomic_alignments ){
	my @algs = $self->fetch_genomic_alignments_for_probe( probe_id=>$rowhash{ probe_id }, chip_type=>$probe->chip_type() );
	$probe->alignments( \@algs );
      }
      push( @probes, $probe );
    }
  }else{
    $self->_pq_fetch_probes_for_gene->execute( $gene_id ) or croak $dbcon->errstr;
    while ( my $hashref = $self->_pq_fetch_probes_for_gene->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      ## doing a little "re-formatting"
      $rowhash{ id } = $rowhash{ probe_id };
      $rowhash{ index_x } = $rowhash{ x };
      $rowhash{ index_y } = $rowhash{ y };
      $rowhash{ adaptor } = $self;
      my %dummyhash = %args;
      $dummyhash{ probe_id } = $rowhash{ probe_id };
      my $probe = CustomCDF::Probe->new( %rowhash );
      if( $load_alignments ){
#	my @algs = $self->fetch_alignments_for_probe_id( probe_id=>$rowhash{ probe_id } );
	my @algs = $self->fetch_alignments_for_probe_id( %dummyhash );
	$probe->alignments( \@algs );
      }
      if( $load_genomic_alignments ){
	my @algs = $self->fetch_genomic_alignments_for_probe( probe_id=>$rowhash{ probe_id }, chip_type=>$probe->chip_type() );
	$probe->alignments( \@algs );
      }
      push( @probes, $probe );
    }
  }
  return( @probes );
}

################################################################
##   Alignment related stuff.
##
################################################################

=head3 fetch_genomic_alignments_for_probe

Fetches genomic alignment(s) for a given probe id.
Note that only "good" alignments are retrieved, i.e. alignments without mismatches!
Internally we have to make shure to just return unique alignments since probe alignments are performed to cDNA and thus we can have multiple alignments of a probe, for more than one transcript of a gene. Converted to chromosomal coordinates, these alignments are then however unique.
Parameters:

=item probe_id      (mandatory) : the probe id.

=item chip_type     (optional)  : the chip_type, either st or at; assumes st if not provided.

=cut
sub fetch_genomic_alignments_for_probe{
   my $self = shift;
   my %args = @_;
   croak "Error: probe_id has to be submitted!" unless exists $args{ probe_id };
   my $chip_type = exists( $args{ chip_type } ) ? $args{ chip_type } : "st";
   my $probe_id = $args{ probe_id };
   ## check first if we do have genomic alignments for the probe.
   ## two ways: 1) easy one: query probe_alignments for seq_type="chromosome",
   ##           2) query probe_to_exon_transcript and calculate genomic coords.
   ## 2) get all alignments to cdna.
   my %AlStrand=();
   my %AlStart=();
   my %AlEnd=();
   my %AlChrom=();
   my $start;
   my $end;
   my $seq_type="chromosome";
   my $dbcon = $self->dbcon;
   $self->_pq_fetch_for_genomic_alignment_for_probe->execute( $probe_id ) or croak $dbcon->errstr;
   while ( my $hashref = $self->_pq_fetch_for_genomic_alignment_for_probe->fetchrow_hashref ){
     my %rowhash = %{$hashref};
     my $alignment_id=$rowhash{ probes_fk };
     ## have to work with hashes in case we do have exon junction probes.
     ## calculate the chromosome_start and end.
     if( "$rowhash{ chromosome_strand }" eq "1" ){
       ## gene is encoded on + strand.
       $start = $rowhash{ exon_chrom_start } + $rowhash{ probe_alignment_exon_start } -1;
       $end = $rowhash{ exon_chrom_start } + $rowhash{ probe_alignment_exon_end } -1;
       if( $chip_type eq "st" ){
	 $AlStrand{ $alignment_id } = "-";
       }else{
	 $AlStrand{ $alignment_id } = "+";
       }
     }else{
       ## gene is encoded on - strand.
       $start = $rowhash{ exon_chrom_end } - $rowhash{ probe_alignment_exon_end } +1;
       $end = $rowhash{ exon_chrom_end } - $rowhash{ probe_alignment_exon_start } +1;
       if( $chip_type eq "st" ){
	 $AlStrand{ $alignment_id } = "+";
       }else{
	 $AlStrand{ $alignment_id } = "-";
       }
     }
     $AlChrom{ $alignment_id } = $rowhash{ chromosome_name };
     if( exists( $AlStart{ $alignment_id } ) ){
       $AlStart{ $alignment_id } = $AlStart{ $alignment_id }.",".$start;
       $AlEnd{ $alignment_id } = $AlEnd{ $alignment_id }.",".$end;
     }else{
       $AlStart{ $alignment_id } = $start;
       $AlEnd{ $alignment_id } = $end;
     }
   }
   ## what we do have now is alignments for unique alignment ids. However, one probe
   ## can have an alignment to each transcript of a gene, which in fact are then a single
   ## genomic alignment. Thus we will reduce this later to single genomic alignments.
   ## 1) get the "real" genomic alignments.
   $self->_pq_fetch_genomic_alignment_for_probe->execute( $probe_id ) or croak $dbcon->errstr;
   while ( my $hashref = $self->_pq_fetch_genomic_alignment_for_probe->fetchrow_hashref ){
     my %rowhash = %{$hashref};
     my $alignment_id=$rowhash{ probes_pk };
     $AlStart{ $alignment_id } = $rowhash{ start };
     $AlEnd{ $alignment_id } = $rowhash{ end };
     $AlStrand{ $alignment_id } = $rowhash{ seq_strand };
     $AlChrom{ $alignment_id } = $rowhash{ seq_name };
   }
   ## ok, and now make the stuff nice and unique.
   my %Alignments = ();
   my $alg_string;
   foreach my $key (keys %AlStart ){
     $alg_string = $AlStart{ $key }.":".$AlEnd{ $key }.":".$AlChrom{ $key }.":".$AlStrand{ $key };
     if( exists( $Alignments{ $alg_string } ) ){
       $Alignments{ $alg_string } = $Alignments{ $alg_string }.",".$key;
     }else{
       $Alignments{ $alg_string } = $key;
     }
   }
   ## ok, and now loop again and make Alignment objects.
   my @splitted;
   my @alignments=();
   foreach my $key (sort( keys %Alignments )){
     @splitted = split( /:/, $key );
     push (@alignments, CustomCDF::Alignment->new( alignment_id=>$Alignments{ $key },
						   probe_id=>$probe_id,
						   start=>[ split( /,/, $splitted[ 0 ] ) ],
						   end=>[ split( /,/, $splitted[ 1 ] ) ],
						   strand=>$splitted[ 3 ],
						   seq_name=>$splitted[ 2 ],
						   seq_type=>"chromosome",
						   mismatches=>0,
						   adaptor=>$self));
   }
   return( @alignments );
}


=head3 fetch_alignments_for_probe_id

Retrieves Alignment(s) for the submitted probe id from the database. Note: only the "probe_alignments" database table is queried unless the option "genomic" is provided.
Parameters and filters:

=item probe_id   (mandatory) : the probe id.

=item seq_name   (optional)  : e.g. seq_name=>"='ENST00000571914'"

=item seq_type   (optional)  : e.g. seq_type=>"='cdna'"

=item start      (optional)  : e.g. start=>"<10000"

=item end        (optional)  : e.g. end=>">=4"

=item seq_strand (optional)  : e.g. seq_strand=>"='+'"

=item mismatches (optional)  : e.g. mismatches=>"=0": defaults to mismatches=>"=0".

=item genomic    (optional)  : fetch all (complete) alignments in genomics coordinates; calls
the fetch_genomic_alignments_for_probe method. This option takes priority over all other options above.

=item chip_type  (optional)  : the chip type, either st or at; see help for fetch_genomic_alignments_for_probe.

=cut
sub fetch_alignments_for_probe_id{
  my $self = shift;
  my %args = @_;
  croak "Error: probe_id has to be submitted!" unless exists $args{ probe_id };
  my $probe_id = $args{ probe_id };
  if( exists( $args{ genomic } ) ){
    return $self->fetch_genomic_alignments_for_probe( %args );
  }
  ## build the query on the fly...
  my $query_string = "select probes_pk, probe_id, seq_name, seq_type, start, end, seq_strand, missmatches from probe_alignments where probe_id='".$probe_id."'";
  ## now checking additional conditions...
  if( exists $args{ seq_name} ){
    $query_string=$query_string." and seq_name".$args{seq_name};
  }
  if( exists $args{ seq_type} ){
    $query_string=$query_string." and seq_type".$args{seq_type};
  }
  if( exists $args{ start} ){
    $query_string=$query_string." and start".$args{start};
  }
  if( exists $args{ end} ){
    $query_string=$query_string." and end".$args{end};
  }
  if( exists $args{ seq_strand} ){
    $query_string=$query_string." and seq_strand".$args{seq_strand};
  }
  if( exists $args{ mismatches} ){
    $query_string=$query_string." and missmatches".$args{mismatches};
  }else{
    $query_string=$query_string." and missmatches=0";
  }
  $query_string=$query_string.";";
  my $db_con=$self->{dbcon};
  my $query = $db_con->prepare( $query_string ) or croak $db_con->errstr;
  $query->execute() or croak $db_con->errstr;
  my @results = ();
  ## looping through the probes and storing them into the array...
  while ( my $hashref = $query->fetchrow_hashref ){
    my %rowhash = %{$hashref};
    $rowhash{alignment_id}=$rowhash{probes_pk};
    $rowhash{mismatches}=$rowhash{missmatches};
    $rowhash{start}=[$rowhash{start}];
    $rowhash{end}=[$rowhash{end}];
    push( @results, CustomCDF::Alignment->new( %rowhash, adaptor=>$self ) );
  }
  # while ( my @row = $query->fetchrow_array){
  #   push( @results, $self->row_to_alignment( @row ) );
  # }
  return( @results );
}


=head1 AUTHOR

Johannes Rainer <johannes.rainer@i-med.ac.at>

=cut

1;

