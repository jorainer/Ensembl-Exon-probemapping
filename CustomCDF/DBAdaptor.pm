package CustomCDF::DBAdaptor;
use strict;
use DBI;
use IO::File;
use warnings;
use Carp;
use Config::Simple;
our $VERSION = "0.3.0";

### package variables:
my $db_con;


=head1 NAME

DBAdaptor - An object to access the database.

=head1 SYNOPSIS

use CustomCDF::DBAdaptor;
my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost", username=>"anonuser", password=>"", dbname=>"homo_sapiens_hugene_cdna_73" );
$db_adaptor::db_con;

=head1 DESCRIPTION

This object provides the connection to a microarray probe alignment and annotation database.

#############
## Database tables and their attributes:
##
# +------------------
# | probe_information
# +------------------
# | - probe_id                    : the id of the probe.
# | - x                           : the x position of the probe on the array.
# | - y                           : the y position of the probe on the array.
# | - platform                    : the microarray platform (microarray name).
# | - sequence                    : the probe's sequence.
# | - nr_chrom_map                : the number of complete genomic alignments.
# | - nr_chrom_map_mm1            : the number of genomic alignments with one mismatch.
# | - nr_chrom_map_mmall          : the number of genomic alignments with one mismatch.
# | - nr_gene_map                 : the number of alignments to different genes.
# | - nr_gene_map_mm              : the number of alignments to genes with a mismatch.
# | - gc_count                    : the number of G-C nucleotides in the probe's sequence.
# | - cel_index                   : the position of the probe within a cel file.
# +------------------
#
# +------------------
# | probe_alignments
# +------------------
# | - probes_pk                   : unique identifier for the alignment.
# | - probe_id                    : the id of the probe.
# | - x                           : the x index of the probe on the array.
# | - y                           : the y index of the probe on the array.
# | - seq_name                    : the name of the sequence, either transcript id or chromosome name.
# | - seq_type                    : the type of sequence, either cdna or chromosome.
# | - start                       : the start position of the alignment.
# | - end                         : the end position of the alignment.
# | - seq_strand                  : the strand of the sequence to which the probe is aligned. Note that for alignments to cdna this should, for sense target arrays always be "-".
# | - missmatches                 : number of mismatches of the alignment.
# +------------------
#
# +------------------
# | probe_to_exon_transcript
# +------------------
# | - transcript_id               : the transcript id.
# | - gene_id                     : the gene id.
# | - probes_fk                   : primary key for the alignment in the probe_alignments table.
# | - chromosome_strand           : the chromosome strand on which the gene is encoded.
# | - probe_id                    : the probe id.
# | - exon_id                     : the exon id.
# | - probe_alignment_exon_start  : start position of the alignment within the exon.
# | - probe_alignment_exon_end    : end position of the alignment within the exon.
# | - position_of_alignment       : for exon junction probes: the order of the alignment segments.
# | - is_junction_alignment       : 1 if the alignment is over an exon junction.
# +------------------

# +------------------
# | exon_transcript
# +------------------
# | - exon_transcript_pk
# | - gene_id
# | - gene_name
# | - gene_biotype
# | - chromosome_name
# | - transcript_id
# | - transcript_chrom_start
# | - transcript_chrom_end
# | - transcript_chrom_strand
# | - transcript_coding_chrom_start
# | - transcript_coding_chrom_end
# | - exon_id
# | - exon_chrom_start
# | - exon_chrom_end
# | - transcript_biotype
# | - coord_system
# +------------------

=head2 Methods

=head3 new

Creates a new connection to the database. host, user, password and db name can be either provided directly or via a configuration file (parameter config_file).

my $db_adaptor = CustomCDF::DBAadaptor->new( host=>"localhost", username=>"anonuser", password=>"", dbname=>"homo_sapiens_hugene_cdna_73" );
$db_adaptor::db_con;

=cut

sub new {
  my($class, %args) = @_;
  my $self = bless({}, $class);
  $self->_init(%args);
  return $self;
}

sub _init{
  my ($self, %args)=@_;
  my $host = exists $args{host} ? $args{host} : '';
  my $username = exists $args{username} ? $args{username} : '';
  my $password = exists $args{password} ? $args{password} : '';
  my $dbname = exists $args{dbname} ? $args{dbname} : '';
  if( exists $args{config_file} ){
    ## load db settings from config file.
    my $cfg = new Config::Simple( $args{config_file} );
    $host=$cfg->param( "database.mysql_host" );
    $username=$cfg->param( "database.mysql_user" );
    $password=$cfg->param( "database.mysql_password" );
    $dbname=$cfg->param( "database.dbprefix" );
  }
  ## connect to DB if everything is there.
  if( $host ne '' & $username ne '' & $dbname ne '' ){
    my $dbcon = DBI->connect( "DBI:mysql:database=$dbname;host=$host:user=$username:password=$password" ) or croak "unable to connect to database $dbname!";
    $self->{dbcon} = $dbcon;
    ## checking database...
    my $table_query = $dbcon->prepare( "show tables" ) or croak $dbcon->errstr;
    $table_query->execute() or croak $dbcon->errstr;
    my $got_probe_information=0;
    my $got_probe_alignments=0;
    my $got_probe_to_exon_transcript=0;
    my $got_exon_transcript=0;
    while ( my @row = $table_query->fetchrow_array ){
      if( $row[ 0 ] eq "probe_information" ){
	$got_probe_information=1;
      }
      if( $row[ 0 ] eq "probe_alignments" ){
	$got_probe_alignments=1;
      }
      if( $row[ 0 ] eq "probe_to_exon_transcript" ){
	$got_probe_to_exon_transcript=1;
      }
      if( $row[ 0 ] eq "exon_transcript" ){
	$got_exon_transcript=1;
      }
    }
    if( ( $got_probe_information + $got_probe_alignments + $got_probe_to_exon_transcript + $got_exon_transcript )!=4 ){
      croak "One or more required database tables are missing!";
    }
  }else{
    croak( "host, username, password and db name have to be provided either directly or using a config_file!" )
  }
  $self->{host} = $host;
  $self->{username} = $username;
  $self->{password} = $password;
  $self->{dbname} = $dbname;
}

=head3 get_probe_adaptor

Returns a CustomCDF::ProbeAdaptor to fetch probes and their alignment from the database. Note: for each ProbeAdaptor a new database connection will be created.

=cut

sub get_probe_adaptor{
  my $self = shift;
  if( !defined $self->{dbcon} ){
    ## this means something went wrong...
    croak( "No database connection available!" )
  }
  my $username = $self->{username};
  my $password = $self->{password};
  my $host = $self->{host};
  my $dbname = $self->{dbname};
  return CustomCDF::ProbeAdaptor->new( host=>$host,
				       username=>$username,
				       password=>$password,
				       dbname=>$dbname
				     );
}

=head3 get_gene_adaptor

Returns a CustomCDF::GeneAdaptor to fetch genes, transcripts and exons from the database. Note: for each GeneAdaptor a new database connection will be created.

=cut

sub get_gene_adaptor{
  my $self = shift;
  if( !defined $self->{dbcon} ){
    ## this means something went wrong...
    croak( "No database connection available!" )
  }
  my $username = $self->{username};
  my $password = $self->{password};
  my $host = $self->{host};
  my $dbname = $self->{dbname};
  return CustomCDF::GeneAdaptor->new( host=>$host,
				       username=>$username,
				       password=>$password,
				       dbname=>$dbname
				     );
}


sub DESTROY{
  ## that's called each time we destroy an instance of Probe.
#  print "byebye\n";
  my $self = shift;
  #print "disconnect from db\n";
  if( defined $self->{dbcon} ){
    $self->{dbcon}->disconnect;
  }
}

=head3 dbcon

returns the db connection.

=cut
sub dbcon{
  my $self=shift;
  unless (ref $self){
    croak( "Can call dbcon only on an object, not a class!" )
  }
  return $self->{dbcon};
}

=head1 AUTHOR

Johannes Rainer <johannes.rainer@i-med.ac.at>

=cut

1;

