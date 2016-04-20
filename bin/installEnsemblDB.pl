#!/usr/bin/perl
## just installing the otherfeatures database.
# arguments:
# -e ensembl version (e.g. 50_36l)
# -d database (e.g. homo_sapiens_otherfeatures_), defaults to homo_sapiens_core_

use DBI;
use Net::FTP;
use Net::SMTP;
use strict;
use Getopt::Std;
use Config::Simple;

my $cfg = new Config::Simple( "../cfg/mysql.cfg" );
### MySQL database settings:
my $mysql_host=$cfg->param( "database.mysql_host" );
my $mysql_username=$cfg->param( "database.mysql_user" );
my $mysql_password=$cfg->param( "database.mysql_password" );
my $mysql_dbname="mysql";
my $mysql_port="3306";

### FTP settings
my $Ensembl_ftp="ftp.ensembl.org";
my $Ensembl_dir="/pub/current_mysql";
my $username="anonymous";
my $password='';   # your email address
###

## where to copy files locally
my $dir="/Volumes/jodata/tmp";

## directory where the scripts are
my $bin_dir="./";

#################################################################################
my $ftp;
my %option=();
#getopts("i:o:f:",\%option);
getopts("e:d:i",\%option);
if( defined( $option{ h } ) ){
  print( "\ninstallEnsemblDB.pl version\n" );
  print( "Downloads and installs a specified Ensembl database locally.\n\n" );
  print( "Usage: perl installEnsemblDB.pl -e:[d:h]\n" );
  print( "Parameters:\n" );
  print( "-e (required): the Ensembl version \n" );
  print( "-d (required): the complete database name, e.g. homo_sapiens_core_82_38. \n" );
  print( "-h : print this help and exit.\n" );
  exit 0;
}

my $version_=$option{ e };
## split the version into the ensembl and database versions.
my @dummy = split /_/,$version_;
my $Ensembl_version=$dummy[0];
my $dbversion=$dummy[1];
$Ensembl_dir="/pub/release-$Ensembl_version/mysql";

##
my $dbname="homo_sapiens_core_";
if( defined( $option{d } )){
    $dbname=$option{ d };
}
installOtherFeatures();

## subroutine for downloading and unzipping the files
sub downloadFiles{
    my @file_array_=@_;
#   print scalar @_."\n\n";
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


##################
## create database
sub createOtherFeatureDB{
  print "creating database...";
  my $dbcon = DBI->connect("dbi:mysql:database=$mysql_dbname;host=$mysql_host;port=$mysql_port","$mysql_username","$mysql_password",{AutoCommit=>0});
  $dbcon->do("CREATE DATABASE $dbname");
  $dbcon->disconnect();
  print "OK\n";
}


sub installOtherFeatures{
  print "installing $dbname database\n";
  $ftp = Net::FTP -> new($Ensembl_ftp);
  $ftp -> login($username,$password) || die "Can't login $!";
  $ftp -> binary;
  $ftp->cwd("$Ensembl_dir/$dbname");
  print "directory on the server: $Ensembl_dir/$dbname\n";
  my @download_files=$ftp->ls();
  my $destdir="$dir/$version_/$dbname";
  ## create the local directory:
  mkdir("$dir/$version_");
  mkdir($destdir);
  downloadFiles(@download_files,$destdir);
  $ftp -> quit();
  ## unzip the files...
  unzipFiles($destdir);
  ## insert them to mySQL
  createOtherFeatureDB();
  ## inserting the database...
  my $theCall = "mysql -u $mysql_username -h $mysql_host --password=$mysql_password $dbname < $destdir/$dbname.sql";
  ## print "Will call: $theCall";
  if(system($theCall)!=0){
    die "error during creation of tables, last message was: $! \n";
  }
  print "importing the data...";
  ## inserting the database...
  if(system("$bin_dir/installDB.sh $destdir $dbname")!=0){
      die "installDB.sh call failed, last message was: $! \n";
    }
  print "finished\n";
  print "Removing database files...";
  if(system("rm --recursive $destdir")!=0){
    die "removing files failed! last message was: $! \n";
  }
  print "finished\n";
}





