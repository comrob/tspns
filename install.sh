#!/bin/sh


# install libraries needed for compiling
sudo apt install -y ccache
sudo apt install libcairo2-dev liblog4cxx-dev libboost-dev

# download and install the algorithm library
sudo apt install git
git clone https://github.com/comrob/crl.git crl
cd crl
./install.sh
cd ..

# download and install GLKH solver
cd glkh
./install.sh
cd ..