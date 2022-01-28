#!/bin/sh

if [ `uname` = FreeBSD ]
then
   mk=gmake
else
   mk=make
fi

# Download and compile LKH Solver

VERSION=1.0

wget -c http://akira.ruc.dk/~keld/research/GLKH/GLKH-1.0.tgz
tar xvfz GLKH-$VERSION.tgz
cd GLKH-$VERSION
$mk
cd -
