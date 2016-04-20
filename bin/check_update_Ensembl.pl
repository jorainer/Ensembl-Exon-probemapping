#!/usr/bin/perl
########################################
##
my $script_version="1.3.0";
# the script checks first if a newer Ensembl version is available online, and if yes it downloads
# and installs the Ensembl_core database. Subsequently it downloads the genomic DNA sequences and
# version 1.1.2 (2011-02-21): downloading "only" the databases and sequences, but NOT starting the alignment.
# version 1.1.3 (2013-07-31): clean up. added options m, a, d, h.
# version 1.1.4 (2013-12-21): downloading also cdna and ncrna fasta files.
# version 1.3.0 (2014-08-07): using git to fetch the perl API instead of CVS (from 76 on).
## usage:
use DBI;
use Net::FTP;
use Net::SMTP;
use strict;
use Getopt::Std;
use Sys::Hostname;
use Config::Simple;

### what do we want to install...
my $install_mysql=0;
my $install_api=0;
my $install_dna=0;

### MySQL database settings:
my $cfg = new Config::Simple( "../cfg/mysql.cfg" );
my $mysql_host=$cfg->param( "mysql_host" );
my $mysql_username=$cfg->param( "mysql_user" );
my $mysql_password=$cfg->param( "mysql_password" );
my $mysql_dbname="mysql";

### SMTP settings
my $send_mail=0;  ## disable the email sending by default.
my $smtp_host;
my $smtp_sender;
my $smtp_to;
############

### FTP settings
my $Ensembl_ftp="ftp.ensembl.org";
my $Ensembl_dir="/pub/current_mysql";
my $username="anonymous";
my $password='johannes.rainer@i-med.ac.at';
###

my $base_dir="/Users/jo/ensembl";
my $unzip=0;

## directory where the scripts are
my $bin_dir="./";

my $species;

### reading parameters.
my %option=();
getopts("admhb:z",\%option);
if( $option{ a } ){
    $install_api=1;
  }
if( $option{ b } ){
  $base_dir=$option{ b };
}
if( $option{ d } ){
    $install_dna=1;
}
if( $option{ m } ){
  $install_mysql=1;
}
if( $option{ z } ){
  $unzip = 1;
}
if( $option{ h } ){
    ## print the help.
    print( "\ncheck_update_Ensembl.pl version ".$script_version."\n" );
    print( "Compares locally installed Ensembl versions with available version online and, if a newer version is available, downloads and installs one or all of API, genomic sequence and MySQL database.\n\n" );
    print( "usage: check_update_Ensembl.pl [-abdhmz]\n" );
    print( "parameters:\n" );
    print( "-h print this help message.\n" );
    print( "-a install the Ensembl perl API.\n" );
    print( "-d download genomic DNA sequence fasta files along with cdna and ncrna fasta files.\n" );
    print( "-m installs the Ensembl core database.\n" );
    print( "-b: provide the path to the local Ensembl installations. Defaults to ~/ensembl.\n\n" );
    print( "-z: automatically unzip the fasta files.\n" );
    exit 0;
}


#################################################################################

## reading folders in the dna_dir, thus we know what versions we have.
opendir(DIR,$base_dir) || die "can not open directory: $base_dir\n";
my @versions = readdir(DIR);
closedir(DIR);

my $ensembl_version_file=$base_dir."/current.cfg";

my $ftp = Net::FTP -> new($Ensembl_ftp);
unless(defined $ftp){
    print "$@\n";
    die "Can't connect!\n";
}
$ftp -> login($username,$password) || die "Can't login $!";
$ftp -> binary;

$ftp->cwd("$Ensembl_dir/");
my $array_ref=$ftp->ls();
my $file;
my $version_;   # combined version string e.g. 43_36e
my $Ensembl_version;
my $NCBI_version;

foreach $file (@$array_ref){
    if($file =~ m/^homo_sapiens_core_(.*)/){
    $version_=$1;
    my @dummy = split /_/,$version_;
    $Ensembl_version=$dummy[0];
    $NCBI_version=$dummy[1];
    }
}
print "check_update_Ensembl.pl version $script_version, Ensembl version on server: $version_ (Ensembl: $Ensembl_version, NCBI: $NCBI_version)\n";

##
my $dna_dir=$base_dir."/".$Ensembl_version."/fasta/";
my $api_dir=$base_dir."/".$Ensembl_version."/API/";

## compare the online version with the versions installed locally:
my $download_=1;
foreach my $dummy (@versions){
    if($dummy eq "." | $dummy eq ".."){
    }
    else{
        print "local version: ".$dummy."\n";
        if($dummy eq $Ensembl_version ){
            $download_=0;
            print "Have already $version_\n";
        }
    }
}



#################### MAIN PART ####################################
## downloading Ensembl core mysql tables
if($download_==1){
  print "Have to fetch new Ensembl version...\n";
  if( $send_mail==1 ){
    ## sending the mail:
    sendMail();
  }
  mkdir( $base_dir."/".$Ensembl_version );

  $species = "homo_sapiens";
  if( $install_mysql==1 ){
    $ftp->cwd("$Ensembl_dir/".$species."_core_".$version_);
    my @download_files=$ftp->ls();
    my $destdir="$base_dir/$Ensembl_version/mysql/".$species."_core_".$version_;
    ## create the local directory:
    mkdir("$base_dir/$Ensembl_version/mysql/");
    mkdir($destdir);
    downloadFiles(@download_files,$destdir);
    $ftp -> quit();
    ## unzip the files...
    unzipFiles($destdir);
    ## insert them to mySQL
    createDB();
    ## inserting the database...
    if(system("$bin_dir/installDB.sh $destdir ".$species."_core_".$version_)!=0){
      die "installDB.sh call failed, last message was: $! \n";
    }
    print "Removing database files...";
    if(system("rm --recursive $destdir")!=0){
      die "removing files failed! last message was: $! \n";
    }
    print "finished\n";
  }

  if( $install_dna==1 | $install_api==1 ){
    if( -e "$base_dir/current" ){
      unlink( "$base_dir/current" );
    }
    system( "ln -s $base_dir/$Ensembl_version $base_dir/current" );
    }

  if( $install_dna==1 ){
    downloadChromosome();
    downloadFasta();
  }
  if( $install_api==1 ){
    ## downloading Ensembl API
    if( $Ensembl_version > 75 ){
      print "Note: fetching the Ensembl Perl API from github\n";
      cloneAPI();
    }else{
      downloadCVS();
    }
  }

  ## download also the mouse core database and FASTA files...
  $species = "mus_musculus";
  if( $install_mysql==1 ){
    $ftp = Net::FTP -> new($Ensembl_ftp);
    $ftp -> login($username,$password) || die "Can't login $!";
    $ftp -> binary;
    $ftp->cwd("$Ensembl_dir/");
    $array_ref=$ftp->ls();
    foreach $file (@$array_ref){
      #    print $file."\n";
      if($file =~ m/^mus_musculus_core_(.*)/){
	$version_=$1;
	my @dummy = split /_/,$version_;
	$Ensembl_version=$dummy[0];
	$NCBI_version=$dummy[1];
	#   print "GOTCHA! ".$file." ".$1."\n";
      }
    }
    $ftp->cwd("$Ensembl_dir/".$species."_core_".$version_);
    my @download_files=$ftp->ls();
    my $destdir="$base_dir/$Ensembl_version/mysql/".$species."_core_".$version_;
    ## create the local directory:
    mkdir($destdir);
    downloadFiles(@download_files,$destdir);
    $ftp -> quit();
    ## unzip the files...
    unzipFiles($destdir);
    createDB();
    ## inserting the database...
    if(system("$bin_dir/installDB.sh $destdir ".$species."_core_".$version_)!=0){
      die "installDB.sh call failed, last message was: $! \n";
    }
    print "Removing database files...";
    if(system("rm --recursive $destdir")!=0){
      die "removing files failed! last message was: $! \n";
    }
    print "finished\n";
  }

  if( $install_dna==1 ){
    downloadChromosome();
    downloadFasta();
  }

  ## write the config file that is used by the mapping scripts!
  open(OUT, "> $ensembl_version_file") or die "cant open data file $ensembl_version_file!";
  print OUT "ensembl=$Ensembl_version\n";
  close(OUT);

}
###############################################################


######################### SUBS ###############################
##
## subroutine for downloading and unzipping the files
sub downloadFiles{
    my @file_array_=@_;
    my $destination_dir_=@_[(scalar @_)-1];
    print "downloading to directory: $destination_dir_\n";
    foreach my $ftp_file (@file_array_){
        if($ftp_file ne $destination_dir_){
            print "downloading $ftp_file...";
            $ftp -> get($ftp_file,"$destination_dir_/$ftp_file");
            print "finished\n";
        }
    }
}

sub unzipFiles{
    my $destination_dir_=@_[0];
    print "unzipping all files...";
    opendir(DIR,$destination_dir_) || die "can not open directory: $destination_dir_\n";
    my @files = readdir(DIR);
    foreach my $file (@files){
        if($file eq "." | $file eq ".." | $file eq "README"){
        }
        else{
            qx(gunzip $destination_dir_"/"$file);
        }
    }
    closedir(DIR);
    print "finished\n";
}


## subroutine to send email
sub sendMail{
    my $smtp = Net::SMTP -> new($smtp_host,Timout=>30,Debug=>0);
    $smtp -> mail($smtp_sender);
    $smtp -> recipient($smtp_to);
    $smtp -> data("Subject: Ensembl update\n\nNew Ensembl version available: $version_\nInstalling on ".hostname.". \nSetting:\n install_api=".$install_api."\n install_dna=".$install_dna."\n install_mysql=".$install_mysql."\n");
    $smtp -> quit();
}

##################
## create database
sub createDB{
    my $dbcon = DBI->connect("dbi:mysql:dbname=$mysql_dbname;host=$mysql_host","$mysql_username","$mysql_password",{AutoCommit=>0});
    $dbcon->do("CREATE DATABASE ".$species."_core_".$version_);
    $dbcon->disconnect();
}

#####
## download and unzip chromosome files.
sub downloadChromosome{
    $ftp = Net::FTP -> new($Ensembl_ftp);
    $ftp -> login($username,$password) || die "Can't login $!";
    $ftp -> binary;
    $ftp->cwd("/pub/current_fasta/$species/dna");
    my @download_files=$ftp->ls();
    my $destdir="$dna_dir/$species/dna/";
    mkdir( $dna_dir );
    mkdir( $dna_dir."/$species" );
    ## create the local directory:
    mkdir($destdir);
    downloadFiles(@download_files,$destdir);
    $ftp -> quit();
    if( $unzip == 1 ){
      unzipFiles($destdir);
    }
}

#####
## download and unzip chromosome files.
sub downloadFasta{
    $ftp = Net::FTP -> new($Ensembl_ftp);
    $ftp -> login($username,$password) || die "Can't login $!";
    $ftp -> binary;
    $ftp->cwd("/pub/current_fasta/$species/cdna");
    my @download_files=$ftp->ls();
    my $destdir="$dna_dir/$species/cdna/";
    mkdir( $dna_dir );
    mkdir( $dna_dir."/$species" );
    ## create the local directory:
    mkdir($destdir);
    downloadFiles(@download_files,$destdir);
    ##$ftp -> quit();
    if( $unzip==1 ){
      ## unzip the files...
      unzipFiles($destdir);
    }
    $ftp->cwd("/pub/current_fasta/$species/ncrna");
    my @download_files=$ftp->ls();
    my $destdir="$dna_dir/$species/ncrna/";
    mkdir( $dna_dir );
    mkdir( $dna_dir."/$species" );
    ## create the local directory:
    mkdir($destdir);
    downloadFiles(@download_files,$destdir);
    $ftp -> quit();
    if( $unzip==1 ){
      ## unzip the files...
      unzipFiles($destdir);
    }
  }


sub downloadCVS{
    ## get Ensembl:
    my $destdir="$api_dir/";
    ## create the local directory:
    mkdir("$api_dir/$Ensembl_version");
    mkdir($destdir);
    print "connecting to ensembl cvs and downloading modules...";
    if(system("$bin_dir/getCVS.sh $destdir $Ensembl_version")!=0){
        die "error downloading ensembl modules! last message was: $! \n";
    }
    print "finished\n";
}

sub cloneAPI{
  my $destdir="$api_dir/";
  ## create the local directory:
  mkdir($destdir);
  ## ensembl core
  if( system( "cd $destdir && git clone https://github.com/Ensembl/ensembl.git" )!=0 ){
    die "error cloning ensembl API! last message was: $! \n";
  }
  ## to be on safe ground... check out the release we want...
  my $thedir = "$destdir/ensembl";
  if( system( "cd $thedir && git checkout release/$Ensembl_version" )!=0 ){
    die "error during release $Ensembl_version check out! Laste message was: $! \n";
  }
  ## ensembl-variation
  if( system( "cd $destdir && git clone https://github.com/Ensembl/ensembl-variation.git" )!=0 ){
    die "error cloning ensembl-variation API! last message was: $! \n";
  }
  ## to be on safe ground... check out the release we want...
  my $thedir = "$destdir/ensembl-variation";
  if( system( "cd $thedir && git checkout release/$Ensembl_version" )!=0 ){
    die "error during release $Ensembl_version check out! Laste message was: $! \n";
  }
  ## ensembl-funcgen
  if( system( "cd $destdir && git clone https://github.com/Ensembl/ensembl-funcgen.git" )!=0 ){
    die "error cloning ensembl-funcgen API! last message was: $! \n";
  }
  ## to be on safe ground... check out the release we want...
  my $thedir = "$destdir/ensembl-funcgen";
  if( system( "cd $thedir && git checkout release/$Ensembl_version" )!=0 ){
    die "error during release $Ensembl_version check out! Laste message was: $! \n";
  }
  ## ensembl-compara
  if( system( "cd $destdir && git clone https://github.com/Ensembl/ensembl-compara.git" )!=0 ){
    die "error cloning ensembl-compara API! last message was: $! \n";
  }
  ## to be on safe ground... check out the release we want...
  my $thedir = "$destdir/ensembl-compara";
  if( system( "cd $thedir && git checkout release/$Ensembl_version" )!=0 ){
    die "error during release $Ensembl_version check out! Laste message was: $! \n";
  }
  print "finished\n";
}



