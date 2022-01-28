#!/bin/sh

cd crl
./clean.sh
cd ..

cd glkh/GLKH-1.0
make clean
cd ../..

make clean 
rm -r include lib