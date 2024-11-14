#!/bin/bash

# Command to run friend generation of hello_world with open data example file
# Run this from the CROWN dir
ANALYSIS=hello_world
CONFIG=config_friends
SAMPLES=dyjets
ERAS=2012
SCOPE=mm
QUANTITIESMAP=${PWD}/build/bin/out_${SCOPE}.root

cd build
# Run cmake to generate c++ code
cmake .. -DANALYSIS=${ANALYSIS} -DCONFIG=${CONFIG} -DSAMPLES=${SAMPLES} -DERAS=${ERAS} -DSCOPES=${SCOPE} -DQUANTITIESMAP=${QUANTITIESMAP}
# Compile generated code
nice -19 make install -j 20
# Use executable on Example file
cd bin
./${CONFIG}_${SAMPLES}_${ERAS}_${SCOPE} out_friend.root ${QUANTITIESMAP}
