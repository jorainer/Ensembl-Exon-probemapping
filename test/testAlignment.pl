use CustomCDF::Alignment;
use strict;
use warnings;

my @starts=( 1, 2, 3 );
my $test = CustomCDF::Alignment->new(
				     start=>\@starts,
				     end=> [ 2, 3, 4 ],
				     strand=>"+",
				     seq_name=>"test",
				     seq_type=>"test"
				    );

print "is gapped? ".$test->is_gapped."\n";
print $test->start;
## $test->start( @starts ); returns error!
$test->start( [ 1, 2, 3 ] );
print "is gapped? ".$test->is_gapped."\n";
print $test->start;

print "the strand: ".$test->strand."\n";
$test->strand( "+" );
print "the strand (expect 1): ".$test->strand."\n";
$test->strand( "-" );
print "the strand (expect -1): ".$test->strand."\n";
$test->strand( -1 );
print "the strand (expect -1): ".$test->strand."\n";
$test->strand( 1 );
print "the strand (expect 1): ".$test->strand."\n";

