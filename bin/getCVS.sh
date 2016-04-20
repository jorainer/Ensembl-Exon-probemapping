#!/bin/sh

DIR=$1
VERSION=$2

mkdir $DIR
cd $DIR
cvs -d :pserver:cvsuser:CVSUSER@cvs.sanger.ac.uk:/cvsroot/CVSmaster login
cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/CVSmaster co -r branch-ensembl-$VERSION ensembl-api

