#!/bin/sh
user=
password=
host=localhost

DIR=$1
DB=$2
DIRDB="$DIR/$DB.sql"

echo "Database: $DB, sqlfile: $DIRDB"

cd $DIR
#mysql -u $user -h $host --password=$password $DB < $DIRDB 1> ../../../log/insert.log 2> ../../../log/insert.err
#mysqlimport -u $user -h $host --password=$password $DB -L *.txt 1> ../../../log/insert.log 2> ../../../log/insert.err
mysql -u $user -h $host --password=$password $DB < $DIRDB
mysqlimport -u $user -h $host --password=$password $DB -L *.txt


