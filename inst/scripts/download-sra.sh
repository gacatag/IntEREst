#!/bin/bash
# Run where your raw files need to be, i.e. RAWDATAPATH
# The SCRIPTSPATH is the path of the scripts folder. 
# For more info see readme.txt. 
cd $RAWDATAPATH
wget -i $SCRIPTSPATH/sraFilesURL.txt

