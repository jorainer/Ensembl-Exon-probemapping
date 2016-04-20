######
my $script_version="0.0.1";
#### Description

use strict;
use Getopt::Std;
use warnings;

my %option=();
getopts( "i:r:h", \%option );
if( defined( $option{ h } ) ){
  print "This is reformat_fasta.pl version ".$script_version."\n\n";
  print "This script reformats a fasta file to make it usable for tools like e.g. RNAplex that have problems with FASTA files in which the sequence is not reported in a single line, but in consecutive lines.\n";
  print "\nParameters:\n";
  print "-i: input file. The output will be written to a file named like the input file, but with a .mod appended to the name.\n";
  print "-r: format sequence id: replace all white spaces in the ID of the sequence with the submitted character. Note, character should be quoted, e.g. \";\"\n";
  print "-h: print this help\n";
}
my $do_format_id = 0;
my $infile;
my $outfile;
my $replace_with;
if( defined( $option{ i } ) ){
  $infile = $option{i};
  my @dummy = split(/\./, $infile );
  my $ending = pop( @dummy );
  push( @dummy, "mod" );
  push( @dummy, $ending );
  $outfile = join( ".", @dummy );
}else{
  die( "input file has to be submitted with -i!" );
}
if( defined( $option{ r } ) ){
  $replace_with = $option{ r };
  $do_format_id = 1;
}


print "reading from file $infile and writing to file $outfile.\n";
if( $do_format_id > 0 ){
  print "Re-formating IDs, replacing white spaces with $replace_with\n";
}

open( IN, "< $infile" ) or die "Can't open input file $infile!\n";
open( OUT, "> $outfile") or die "Can't open output file $outfile for writing!\n";
my $sequence="empty";
my $id = "empty";
while( <IN> ){
  chomp;
  if( /^>/ ){
    if( $sequence ne "empty" ){
      print OUT $sequence."\n";
    }
    $id = $_;
    if( $do_format_id > 0 ){
      $id =~ s/\s/$replace_with/g;
    }
    print OUT $id."\n";
    $sequence="";
  }else{
    $sequence=$sequence.$_;
  }
}
## write the last one...
print OUT $sequence."\n";

close( IN );
close( OUT );

