VCFs to SNP alignment
=====================

USDA APHIS Veterinary Services (VS) Mycobacterium tuberculosis complex, mainly M. bovis, and Brucella sp. genotyping from whole genome sequence (WGS) outputting BAM, VCF, SNP tables and phylogentic trees. 

SUMMARY
========

Using BWA alignments VCFs are created for each isolate using GATK's haplotype caller.  The overview below discribes a two step pipeline, briefly: step 1 - outputing VCFs and step 2 - gathering those VCFs to produce SNP alignments, tables and phylogenetic trees of grouped isolates.  This two step pipeline has been developed to provide routine diagnostics results rapidly, with quality assurance and easy to intrepret reporting.

## Script 1
Step 1 is fairly straight forward.  Our main workflow includes Mycobacterium tuberculosis complex and Brucella sp. and therefore the script has been optimized for such.  The script begins by selecting the best reference and determining the spoliogtype for TB complex isolates and MLST for Brucella spieces.  The script then goes on to using what is now a fairly standard bacterial SNP calling pipeline where BWA aligns and GATK calls SNPs.  However, in addition to what is output by the HaplotypeCaller, Map zero positions are added to the VCFs.  This is not a standard GATK option and is done using shell commands.  Including Map zero positions allows a more accurate SNP summary to be represented in step 2.

Script 1 is at:  https://github.com/USDA-VS/snp_analysis/blob/master/bin/processzips.sh

## Script 2
Step 2 is called on VCFs output from step 1.  References chosen in script 1 have been selected because they have been found to be relatively close to the isolate.  The closer the reference is to the isolate the less overall SNP calling error is seen.  VCFs analyzed in step 2 must all be output from the same reference.  Obviously VCFs analyzed using different references can not be used in the same comparison.

In addition to choosing a closely related reference, which minimizes SNP calling error, there are three additional external files, or dependencies, used to create high quality, informative SNP alignments.  As shown in bovis_dependency_view.jpg and suis1_dependency_view.jpg the reference used to build VCFs must be reflected in the three dependencies.  The three dependent files are: filter file, defining SNPs, and gbk file.

INSTALL - Script2
=================

## Python environment setup

Instructions shown are to install with user, not root, privileges.  Use root privileges if available or with another user PATH setup if desired.

Script 2 is written in Python and must be ran using Python 3.  

Anaconda is a highly trusted Python package distrubution platform.  If running Python 2 a new environment can be set without disrupting your current Python environment.  See note below for installing an additional Anaconda environment.  

Install Anaconda if not already installed.  Tested using Anaconda3-4.3.1, or try the latest at: https://www.continuum.io/downloads

On Mac OS X

    ~$ wget https://repo.continuum.io/archive/Anaconda3-4.3.1-MacOSX-x86_64.sh
    ~$ bash Anaconda3-4.3.1-MacOSX-x86_64.sh

On Linux

    ~$ wget https://repo.continuum.io/archive/Anaconda3-4.3.1-Linux-x86_64.sh        
    ~$ bash Anaconda3-4.3.1-Linux-x86_64.sh
    
Use Anaconda's default installation except when asked if to, "prepend to PATH", choose yes.
    
Once Anaconda is installed setup Bioconda channels.  Add them in the order shown below.  Order is important.

    ~$ conda config --add channels conda-forge
    ~$ conda config --add channels defaults
    ~$ conda config --add channels r
    ~$ conda config --add channels bioconda
    
As of Anaconda3-4.3.1 ete3 requires python version < 3.6

    ~$ conda install python=3.5
    
    ~$ conda install ete3 pyvcf biopython
    ~$ conda update ete3 pyvcf biopython

Xvfb (short for X virtual framebuffer) is a display server implementing the X11 display server protocol and must be in your environment to generate PDF and SVG tree files.  If X11 is not in your environment xvfbwrapper will install but not work.  Xvfb is likely already install but beware.  Root privileges will be needed if it is not yet in your environment.  It's available on Mac OS X via XQuartz and Linux via your package manager.

Note, Xvfb is finicky.  Using Xvfb is an unfortuante necessity.  Any feedback to improve the portablity generating PDF and SVG files will be appreciated.   On Mac OS X install using pip and on Linux use conda.

On Mac OS X

Install XQuartz https://www.xquartz.org then:

    ~$ pip install xvfbwrapper

On Linux

    ~$ sudo apt-get install xvfb
    ~$ conda install xvfbwrapper
    
pandas 0.18.1 is required as of pandas 0.20.0

    ~$ conda install pandas=0.18.1

RAxML must be in your PATH as: raxmlHPC-SSE3, raxmlHPC-PTHREADS-AVX2, or raxml.  It seems raxmlHPC-SSE3 is the most system universal but if you have the correct computer architecture running raxmlHPC-PTHREADS-AVX2 is faster.  The script will first look for raxmlHPC-PTHREADS-AVX2.  If it is not found it will look for raxmlHPC-SSE3, then raxml.  If none are found in your PATH the script will fail.

If RAxML is not already installed the fastest way to install is:

    ~$ conda install raxml
    ~$ conda update raxml

## Script and file dependents
The script will look in your home directory for file dependencies.  

    $ ~/dependencies
    
To install dependencies, run the command below with your home directory set as your current working directory.  Check this repo periodically for updates.

Clone dependencies

    ~$ git clone https://github.com/stuber/dependencies.git

Clone script: 

    ~$ git clone https://github.com/stuber/VCFs_to_SNP_alignment.git

Change directory to `VCFs_to_SNP_alignment` and run line below to put script in your PATH.

    $ thepath=$(pwd); ln -s ${thepath}/script2.py ~/anaconda3/bin
    
## Test
Use files bundled with dependencies to test.  Make working directory that containing VCFs and call script.  With the command above your script should be in your PATH.

For list of options:
    
    $ script2.py -h
    
If running bovis VCFs, run the following:

    $ script2.py -s bovis

Debug with

    $ script2.py -s bovis -d

Adding an additional environment
=======================================

If a new Anaconda environment is needed without making changes to your current:
        
    $ conda create -n anaconda400 anaconda=4.0.0 anaconda
    
To activate this environment, use:
    
    > source activate anaconda400
    
To deactivate this environment, use:
    
    > source deactivate anaconda400

OVERVIEW
========

*WORK IN PROGRESS*

Importance of step 1.  This is where the SNPs are called.  The final analysis is only going to be as good as the SNPs being called.  It doesn't matter which programs are used, per se, as long as each SNP called (or not called) can be justified and validated.

Because of this step 1 and step 2 are very much linked to each other.


OBJECTIVE
==========

