package CustomCDF::GeneAdaptor;
use strict;
use DBI;
use IO::File;
use warnings;
use Carp;
use Config::Simple;
use CustomCDF::Feature;
use CustomCDF::Exon;
use CustomCDF::Transcript;
use CustomCDF::Gene;
our $VERSION = "0.3.0";

our @ISA = qw( CustomCDF::DBAdaptor );

my $query_fetch_gene_by_id="select distinct gene_id, gene_name, gene_biotype, chromosome_name, transcript_chrom_strand, coord_system from exon_transcript where gene_id=?;";

my $query_fetch_exon_by_id="select exon_id, exon_chrom_start, exon_chrom_end, chromosome_name, transcript_chrom_strand, coord_system from exon_transcript where exon_id=?;";

my $query_fetch_transcript="select distinct transcript_id, transcript_biotype, transcript_chrom_start, transcript_chrom_end, transcript_chrom_strand, transcript_coding_chrom_start, transcript_coding_chrom_end, chromosome_name, coord_system from exon_transcript where transcript_id=?";

my $query_fetch_transcript_with_exon="select distinct transcript_id, transcript_biotype, transcript_chrom_start, transcript_chrom_end, transcript_chrom_strand, transcript_coding_chrom_start, transcript_coding_chrom_end, chromosome_name, coord_system, exon_id, exon_chrom_start, exon_chrom_end from exon_transcript where transcript_id=? order by exon_chrom_start";

my $query_fetch_gene_by_transcript="select distinct gene_id, gene_name, gene_biotype, chromosome_name, transcript_chrom_strand, coord_system from exon_transcript where transcript_id=?;";

my $query_fetch_gene_with_transcript="select distinct gene_id, gene_name, gene_biotype, chromosome_name, transcript_id, transcript_chrom_start, transcript_chrom_end, transcript_chrom_strand, transcript_coding_chrom_start, transcript_coding_chrom_end, transcript_biotype, coord_system from exon_transcript where gene_id=? order by transcript_chrom_start;";

my $query_fetch_exons_for_gene="select distinct exon_id, exon_chrom_start, exon_chrom_end, transcript_chrom_strand, chromosome_name, coord_system from exon_transcript where gene_id=? order by exon_chrom_start";

my $query_fetch_exons_for_transcript="select distinct exon_id, exon_chrom_start, exon_chrom_end, transcript_chrom_strand, chromosome_name, coord_system from exon_transcript where transcript_id=? order by exon_chrom_start";


=head1 NAME

GeneAdaptor - An object to access genes stored in the database.

=head1 SYNOPSIS

use CustomCDF::DBAdaptor;
my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost", username=>"anonuser", password=>"", dbname=>"homo_sapiens_hugene_cdna_74" );
my $gene_adaptor = $db_adaptor->get_gene_adaptor();

=head1 DESCRIPTION

This object allows to fetch genes along with their transcripts and exons from the alignment database.
CustomCDF::GeneAdaptor inherits from CustomCDF::DBAdaptor directly.

=head2 Methods

=head3 new

=cut
sub new{
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
  my $pq_fetch_gene_by_id=$dbcon->prepare_cached( $query_fetch_gene_by_id ) or croak $dbcon->errstr;
  $self->{pq_fetch_gene_by_id}=$pq_fetch_gene_by_id;
  my $pq_fetch_exon_by_id=$dbcon->prepare_cached( $query_fetch_exon_by_id ) or croak $dbcon->errstr;
  $self->{pq_fetch_exon_by_id}=$pq_fetch_exon_by_id;
  my $pq_fetch_gene_with_transcript=$dbcon->prepare_cached( $query_fetch_gene_with_transcript ) or croak $dbcon->errstr;
  $self->{pq_fetch_gene_with_transcript}=$pq_fetch_gene_with_transcript;
  my $pq_fetch_gene_by_transcript=$dbcon->prepare_cached( $query_fetch_gene_by_transcript ) or croak $dbcon->errstr;
  $self->{pq_fetch_gene_by_transcript}=$pq_fetch_gene_by_transcript;
  my $pq_fetch_transcript=$dbcon->prepare_cached( $query_fetch_transcript ) or croak $dbcon->errstr;
  $self->{pq_fetch_transcript}=$pq_fetch_transcript;
  my $pq_fetch_transcript_with_exon=$dbcon->prepare_cached( $query_fetch_transcript_with_exon ) or croak $dbcon->errstr;
  $self->{pq_fetch_transcript_with_exon}=$pq_fetch_transcript_with_exon;
  my $pq_fetch_exons_for_gene=$dbcon->prepare_cached( $query_fetch_exons_for_gene ) or croak $dbcon->errstr;
  $self->{pq_fetch_exons_for_gene}=$pq_fetch_exons_for_gene;
  my $pq_fetch_exons_for_transcript=$dbcon->prepare_cached( $query_fetch_exons_for_transcript ) or croak $dbcon->errstr;
  $self->{pq_fetch_exons_for_transcript}=$pq_fetch_exons_for_transcript;
}

#####
## getters for the cached queries.
sub _pq_fetch_gene_by_id{
  my $self = shift;
  return( $self->{pq_fetch_gene_by_id} );
}
sub _pq_fetch_exon_by_id{
  my $self = shift;
  return( $self->{pq_fetch_exon_by_id} );
}
sub _pq_fetch_transcript{
  my $self = shift;
  return( $self->{pq_fetch_transcript} );
}
sub _pq_fetch_transcript_with_exon{
  my $self = shift;
  return( $self->{pq_fetch_transcript_with_exon} );
}
sub _pq_fetch_gene_with_transcript{
 my $self = shift;
 return( $self->{pq_fetch_gene_with_transcript} );
}
sub _pq_fetch_gene_by_transcript{
 my $self = shift;
 return( $self->{pq_fetch_gene_by_transcript} );
}
sub _pq_fetch_exons_for_gene{
 my $self = shift;
 return( $self->{pq_fetch_exons_for_gene} );
}
sub _pq_fetch_exons_for_transcript{
 my $self = shift;
 return( $self->{pq_fetch_exons_for_transcript} );
}


sub DESTROY{
  ## that's called each time we destroy an instance of Probe.
#  print "byebye\n";
  my $self = shift;
#  $self->{dbcon}->disconnect;
  ## finishing all
  #print "finishing cached queries\n";
  if( defined $self->dbcon ){
    $self->_pq_fetch_gene_by_id->finish;
    $self->_pq_fetch_exon_by_id->finish;
    $self->_pq_fetch_gene_with_transcript->finish;
    $self->_pq_fetch_gene_by_transcript->finish;
    $self->_pq_fetch_exons_for_gene->finish;
    $self->_pq_fetch_transcript->finish;
    $self->_pq_fetch_exons_for_transcript->finish;
    $self->_pq_fetch_transcript_with_exon->finish;
  }
}


################################################################
##   Gene related stuff.
##
################################################################

=head3 fetch_gene

That's supposed to be the super-duper method that allows to fetch a gene with all its transcripts and all the transcripts exons.
Parameters:

=item gene_id          (mandatory): the gene id.

=item load_transcripts (optional) : whether transcripts of the gene should also be fetched.

=item load_exons       (optional) : whether transcripts and exons should also be fetched.

Example:
use CustomCDF::GeneAdaptor;
use CustomCDF::DBAdaptor;
use CustomCDF::Gene;

my $db_adaptor = CustomCDF::DBAdaptor->new( host=>"localhost",
						  username=>"anonuser",
						  password=>"",
						  dbname=>"homo_sapiens_hugene_cdna_74"
				     );
my $gene_adaptor = $db_adaptor->get_gene_adaptor();

my $gene = $gene_adaptor->fetch_gene( gene_id=>"ENSG00000000003" );


=cut
sub fetch_gene{
  my $self = shift;
  ## here we assume some additional stuff...
  my %args = @_;
  croak "Error: the gene id has to be speficied with gene_id=>!" unless exists( $args{ gene_id } );
  my $gene_id = $args{gene_id};
  my $load_transcripts = exists( $args{ load_transcripts } );
  my $load_exons = exists( $args{ load_exons } );
  my $dbcon = $self->dbcon;
  ## do different stuff depending on the setup...
  ## plainly query gene and return that.
  if( !$load_transcripts & !$load_exons ){
    $self->_pq_fetch_gene_by_id->execute($gene_id) or croak $dbcon->errstr;
    my %rowhash = %{ $self->_pq_fetch_gene_by_id->fetchrow_hashref };
    if( %rowhash ){
      return CustomCDF::Gene->new( %rowhash, adaptor=>$self );
    }
  }
  ## load also the transcripts.
  if( $load_transcripts & !$load_exons ){
    $self->_pq_fetch_gene_with_transcript->execute($gene_id) or croak $dbcon->errstr;
    my $do_gene=1;
    my $do_have_result=0;
    my $gene;
    my @transcripts = ();
    while ( my $hashref = $self->_pq_fetch_gene_with_transcript->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      if( $do_gene==1 ){
    	$do_gene=0;
    	$do_have_result=1;
    	$gene = CustomCDF::Gene->new( %rowhash, adaptor=>$self );
      }
      ## and the transcripts...
      push( @transcripts, CustomCDF::Transcript->new( %rowhash, adaptor=>$self ) );
    }
    if( $do_have_result==1 ){
      ## add the transcripts...
      $gene->transcripts( \@transcripts );
      return( $gene );
    }
  }
  if( $load_exons ){
    ## do the full stuff. also do load transcripts.
    ## we query distinct gene and transcript ids from the exon_transcript.
    ## for the first row we define the gene.
    ## for each row we call fetch_transcript. we have thus more queries.
    $self->_pq_fetch_gene_with_transcript->execute($gene_id) or croak $dbcon->errstr;
    my $do_gene=1;
    my $do_have_result=0;
    my $gene;
    my @transcripts = ();
    while ( my $hashref = $self->_pq_fetch_gene_with_transcript->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      if( $do_gene==1 ){
    	$do_gene=0;
    	$do_have_result=1;
    	$gene = CustomCDF::Gene->new( %rowhash, adaptor=>$self );
      }
      ## and the transcripts: query again!
      push( @transcripts, $self->fetch_transcript( transcript_id=>$rowhash{transcript_id}, load_exons=>1 ) );
    }
    if( $do_have_result==1 ){
      ## add the transcripts...
      $gene->transcripts( \@transcripts );
      return( $gene );
    }
  }
}


=head3 fetch_transcripts_for_gene

Loads all transcripts for a given gene.
Parameters:

=item gene_id     (mandatory)  : the gene id.

=item load_exons  (optional)   : load also the exons for each transcript.

=cut
sub fetch_transcripts_for_gene{
  my $self = shift;
  ## here we assume some additional stuff...
  my %args = @_;
  croak "Error: the gene id has to be speficied with gene_id=>!" unless exists( $args{ gene_id } );
  $args{ load_transcripts }=1;
  my $gene = $self->fetch_gene( %args );
  return $gene->transcripts();
}

=head3 fetch_exons_for_gene

Loads all (unique) Exons for a given gene from the database.
The difference to "load_exons" parameter in the fetch_gene is that we are here fetching all exons without their transcripts! Thus we do have an array with unique Exons, ordered by chromosome start.

=cut
sub fetch_exons_for_gene{
  my $self = shift;
  ## here we assume some additional stuff...
  my %args = @_;
  croak "Error: the gene id has to be speficied with gene_id=>!" unless exists( $args{ gene_id } );
  my $gene_id = $args{gene_id};
  my $dbcon = $self->dbcon;
  $self->_pq_fetch_exons_for_gene->execute( $gene_id ) or croak $dbcon->errstr;
  my @exons = ();
  while ( my $hashref = $self->_pq_fetch_exons_for_gene->fetchrow_hashref ){
    my %rowhash = %{$hashref};
    push( @exons, CustomCDF::Exon->new( %rowhash, adaptor=>$self ) );
  }
  return @exons;
}


################################################################
##   Transcript related stuff.
##
################################################################

=head3 fetch_transcript

Load a Transcript object from the database for a given transcript id.
Parameters:

=item transcript_id (mandatory) : the id of the transcript.

=item load_exons    (optional)  : load also all exons for the gene. Exons will be retrieved ordered by their chromosomal start position.

=cut
sub fetch_transcript{
  ## two possibilities, just the transcript, or the transcript with the exons.
  my $self = shift;
  my %args = @_;
  croak "Error: the transcript id has to be speficied with transcript_id=>!" unless exists( $args{ transcript_id } );
  my $transcript_id = $args{transcript_id};
  my $load_exons = exists( $args{ load_exons } );
  my $dbcon = $self->dbcon;
  if( !$load_exons ){
    ## return just the transcript
    $self->_pq_fetch_transcript->execute($transcript_id) or croak $dbcon->errstr;
    my %rowhash = %{ $self->_pq_fetch_transcript->fetchrow_hashref };
    if( %rowhash ){
      my $transcript = CustomCDF::Transcript->new( %rowhash, adaptor=>$self );
      return $transcript;
    }
  }
  if( $load_exons ){
    $self->_pq_fetch_transcript_with_exon->execute($transcript_id) or croak $dbcon->errstr;
    my $do_transcript=1;
    my $do_have_result=0;
    my $transcript;
    my @exons = ();
    while ( my $hashref = $self->_pq_fetch_transcript_with_exon->fetchrow_hashref ){
      my %rowhash = %{$hashref};
      if( $do_transcript==1 ){
    	$do_transcript=0;
    	$do_have_result=1;
    	$transcript = CustomCDF::Transcript->new( %rowhash, adaptor=>$self );
      }
      ## and the transcripts...
      push( @exons, CustomCDF::Exon->new( %rowhash, adaptor=>$self ) );
    }
    if( $do_have_result==1 ){
      ## add the transcripts...
      $transcript->exons( \@exons );
      return( $transcript );
    }
  }
}


=head3 fetch_exons_for_transcript

Loads all Exons for a given transcript

=cut
sub fetch_exons_for_transcript{
  my $self = shift;
  ## here we assume some additional stuff...
  my %args = @_;
  croak "Error: the transcript id has to be speficied with transcript_id=>!" unless exists( $args{ transcript_id } );
  my $transcript_id = $args{transcript_id};
  my $dbcon = $self->dbcon;
  $self->_pq_fetch_exons_for_transcript->execute( $transcript_id ) or croak $dbcon->errstr;
  my @exons = ();
  while ( my $hashref = $self->_pq_fetch_exons_for_transcript->fetchrow_hashref ){
    my %rowhash = %{$hashref};
    push( @exons, CustomCDF::Exon->new( %rowhash, adaptor=>$self ) );
  }
  return @exons;
}


=head3 fetch_gene_id_by_transcript

Retrieves the gene id for a given transcript id.

=cut
sub fetch_gene_id_by_transcript{
  my $self = shift;
  ## here we assume some additional stuff...
  my %args = @_;
  croak "Error: the transcript id has to be speficied with transcript_id=>!" unless exists( $args{ transcript_id } );
  my $transcript_id = $args{transcript_id};
  my $dbcon = $self->dbcon;
  $self->_pq_fetch_gene_by_transcript->execute( $transcript_id ) or croak $dbcon->errstr;
  my @results = ();
  while ( my $hashref = $self->_pq_fetch_gene_by_transcript->fetchrow_hashref ){
    my %rowhash = %{$hashref};
    push( @results, $rowhash{ gene_id } );
  }
  if( scalar( @results )==1 ){
    return $results[ 0 ];
  }
  return @results;
}

=head3 fetch_gene_by_transcript

Retrieves the gene for a given transcript id.

=cut
sub fetch_gene_by_transcript{
  my $self = shift;
  ## here we assume some additional stuff...
  my %args = @_;
  croak "Error: the transcript id has to be speficied with transcript_id=>!" unless exists( $args{ transcript_id } );
  my $transcript_id = $args{transcript_id};
  my $dbcon = $self->dbcon;
  $self->_pq_fetch_gene_by_transcript->execute( $transcript_id ) or croak $dbcon->errstr;
  my @results = ();
  while ( my $hashref = $self->_pq_fetch_gene_by_transcript->fetchrow_hashref ){
    my %rowhash = %{$hashref};
    push( @results, CustomCDF::Gene->new( %rowhash, adaptor=>$self ) );
  }
  if( scalar( @results )==1 ){
    return $results[ 0 ];
  }
  return @results;
}


################################################################
##   Exon related stuff.
##
################################################################

=head3 fetch_exon

Load an Exon object from the database for a given transcript id.
Parameters:

=item exon_id   (mandatory)   : the id of the exon.

=cut
sub fetch_exon{
  my $self = shift;
  my %args = @_;
  croak "Error: the exon id has to be speficied with exon_id=>!" unless exists( $args{ exon_id } );
  my $id = $args{exon_id};
  my $dbcon = $self->dbcon;
  $self->_pq_fetch_exon_by_id->execute($id) or croak $dbcon->errstr;
  my %rowhash = %{ $self->_pq_fetch_exon_by_id->fetchrow_hashref };
  if( %rowhash ){
    my $exon = CustomCDF::Exon->new( %rowhash, adaptor=>$self );
    return $exon;
  }
}


1;
