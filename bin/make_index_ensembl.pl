#######
my $script_version_="0.2.1";
##
# script to build bowtie indices from genomic FASTA files downloaded from Ensembl
## options:
# -d directory where the Ensembl fasta files can be found.
#######
use IO::File;
use Getopt::Std;
use strict;
use warnings;

### settings
my $bowtie_bin="bowtie-build";
my $bowtie_options="--seed 123";
my $bowtie2_bin="bowtie2-build";
my $bowtie2_options="--seed 123 --large-index";
my $gmap_bin="gmap_build";
my $gmap_options="";
## other stuff
my $directory=".";
my $species="Homo_sapiens";
my $aligner="bowtie";
my $file_format=2;
my $ensembl_release="unknown";

my %option=();
getopts("a:d:e:f:s:h",\%option);

### print the help message:
if( defined( $option{h} ) ){
  print( "\nmake_index_ensembl.pl version ".$script_version_.".\n" );
  print( "Creates an index required for fast sequence alignment using the bowtie or gmap/gsnap aligner.\n\n" );
  print( "usage: make_index_ensembl.pl.pl -[ad:hs:]\n" );
  print( "parameters:\n" );
  print( "-a optional: the aligner for which the index should be generated.\n" );
  print("-e optional: the Ensembl release. This will be used in the generated output file name. Note: this is mandatory!\n");
  print( "-d optional: the path to the fasta files with the genomic sequence. By default, the script looks for fasta files in the current directory.\n" );
  print( "-f optional: define the file name format of the genome fasta files. Options are 1 or 2 (the default), Ensembl releases up to Ensembl 75 (option 1) included the Ensembl version name in the file name, after release 75 (option 2) the version name is no longer included.\n" );
  print( "-h : print this help message.\n");
  print( "-s optional: specify the species; defaults to Homo_sapiens. Note that the species should match the species name given in the file names of the fasta file with the genomic sequences.\n\n");
  exit 0;
}

if(defined($option{e})){
  $ensembl_release = $option{e};
}
if(defined($option{f})){
  my $fformat = $option{f};
  if($fformat eq "1"){
    $file_format = 1;
  }elsif($fformat eq "2"){
    $file_format = 2;
  }elsif($fformat == 1){
    $file_format = 1;
  }elsif($fformat == 2){
    $file_format = 2;
  }
  if($file_format!=1 & $file_format!=2){
    die("Only 1 or 2 are allowed as file format values!");
  }
}
if( defined( $option{d} ) ){
	$directory = $option{d};
}
if( defined( $option{s} ) ){
	$species = $option{s};
}
if( defined( $option{a} ) ){
  $aligner = $option{a};
  if( $aligner ne "bowtie" ){
    if( $aligner ne "gmap" ){
      if( $aligner ne "bowtie2" ){
	stop( "only bowtie, bowtie2 or gmap are allowed for parameter -a!" );
      }
    }
  }
}


print "this is make_bowtie_index_ensembl.pl version ".$script_version_."\n";

## todo:
# read all files in directory,
# loop through files, if ok add to a "," list
# start bowtie_build

## reading the files in this directory:
opendir(DIR,$directory) || die "can not open directory: $directory\n";
my @all_files = readdir(DIR);
closedir(DIR);

my $assembly;
my $release;
my @fastafiles=();

## assuming Ensembl format: <species>.<assembly>.<release>.<sequence type>.<id type>.<id>.fa
## e.g. Homo_sapiens.GRCh37.58.dna.chromosome.13.fa
my $idx_assembly = 1;
my $idx_release = 2;  ## Ensembl release... only for Ensembl pre 76
my $idx_type = 3;
my $idx_seqtype = 4;
my $idx_chromno = 5;
if($file_format==2){
  ## the file format without the Ensembl release:
  ## e.g. Homo_sapiens.GRCh38.dna.chromosome.13.fa
  $idx_type=2;
  $idx_seqtype=3;
  $idx_chromno=4;
}
foreach my $file ( @all_files ){
#    print $file."\n";
    if($file =~ m/^$species(.*)/){
      my @dummy = split /\./,$file;
      $assembly = $dummy[$idx_assembly];
      $release = $dummy[$idx_release];
      ## check if the file is a chromosome fasta file (no repeat masked (rm))
      if( $dummy[$idx_type] eq "dna" & $dummy[$idx_seqtype] eq "chromosome" ){
	## have to restrict to the default "human chromosomes"
	if( $dummy[$idx_chromno] =~ m/^(\d*|X|Y|MT)$/ ){
	  # if( $filestring eq "empty" ){
	  #   $filestring = $directory."/".$file;
	  # }
	  # else{
	  #   $filestring = $filestring.",".$directory."/".$file;
	  # }
	  push( @fastafiles, $directory."/".$file );
	  ## check if that file is gzipped:
	  if( $file =~ m/(.*).gz$/ ){
	    $gmap_options="-g";
	  }
	}
      }
    }
}

if( scalar( @fastafiles )==0 ){
  die "no files in the correct fasta file format found in directory ".$directory."\nIf the fasta files are from Ensembl versions prior 76 you should add -f 1.\n";
}
if($file_format==1){
  $ensembl_release = $release;
}

if( $aligner eq "bowtie" ){
  my $output = $species."_".$assembly."_".$ensembl_release;
  my $filestring = join( ",", sort( @fastafiles ) );
  print "running ".$bowtie_bin." ".$bowtie_options." ".$filestring." ".$output."\n";
  if( system( $bowtie_bin." ".$bowtie_options." ".$filestring." ".$output )!=0 ){
    die "ERROR, bowtie call failed!\n";
  }
}
if( $aligner eq "bowtie2" ){
  my $output = $species."_".$assembly."_".$ensembl_release;
  my $filestring = join( ",", sort( @fastafiles ) );
  print "running ".$bowtie2_bin." ".$bowtie2_options." ".$filestring." ".$output."\n";
  if( system( $bowtie2_bin." ".$bowtie2_options." ".$filestring." ".$output )!=0 ){
    die "ERROR, bowtie2 call failed!\n";
  }
}
if( $aligner eq "gmap" ){
  my $output = "-d ".$species."_".$assembly."_".$ensembl_release;
  my $filestring = join( " ", sort( @fastafiles ) );
  print "running ".$gmap_bin." ".$gmap_options." ".$output." ".$filestring."\n";
  if( system( $gmap_bin." ".$gmap_options." ".$output." ".$filestring )!=0 ){
    die "ERROR, gmap call failed!\n";
  }
}


