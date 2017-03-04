This README.TXT file is located in the scripts/ folder of the
MDS.Chr22 package. This folder contains the scripts that
were used to produce the data located in the extdata/ folder.

This file describes how to regenrate the data in extdata/
folder.

-----
Requirements:

- The samples were sequenced very deep (~ 60-80 million pair of 
sequenced reads) hence the mapping could be extremely time-consuming (from 
24 hours up to 5 days). In this respect, we recommend running the mapping 
(Tophat2 runs) of the sequencing reads by using at least 4 threads and in 
an environment with at least 4 computing cores. We ran the analysis
on Introni i.e. a local server with 64 cores (AMD Opteron(TM) Processor 6274)
528 Gb of RAM, and a Linux x86_64 (CentOS release 5.11) Operating System.

- 1.3 TB of available disk space.

- Fast internet connection. About 150 GB needs to be downloaded. 

- Software: bowtie2-2.0.0-beta7, tophat-2.0.5, and SAMtools (0.1.18).

-----
Installing required software:

- For Bowtie2 instllation see 
   http://bowtie-bio.sourceforge.net/tutorial.shtml

- For Tophat2 installation see 
   https://ccb.jhu.edu/software/tophat/tutorial.shtml

- For samtools installation see
   http://www.htslib.org/download/

-----
Regenrating data:
To regenerate the data carry out the following: 

1. Set the proper environmental variables: SCRIPTSPATH should be set to the 
   absolute path of scripts/ folder of IntEREst package. Set RAWDATAPATH to 
   where you like to download the raw files. If IntEREst is installed, in R you 
   can retreive the path to the scripts with the following command:
   system.file("scripts", package="IntEREst"). Set MAPPINGPATH to where the 
   mapping results and the contents in extdata/ should be stored (Could be the 
   same as RAWDATAPATH). indexPATH should be set as the path of the bowtie 2 
   indexes.
    

2. Run the commands in following order:

   sh $SCRIPTSPATH/download-sra.sh
   sh $SCRIPTSPATH/sraTofastqConvert.sh
   sh $SCRIPTSPATH/allTophat2Runs.sh
   sh $SCRIPTSPATH/selectCh22fromBam.sh
   Rscript --vanilla $SCRIPTSPATH/build_mdsChr22Obj.R $MAPPINGPATH ./

