use DBI;
use Net::FTP;
use Net::SMTP;
use strict;
use Getopt::Std;
use Config::Simple;
my $dbcon = DBI->connect("dbi:mysql:database=mysql;host=localhost;port=3306","jo","jo123",{AutoCommit=>0});


