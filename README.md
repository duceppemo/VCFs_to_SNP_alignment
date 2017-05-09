VCFs to SNP alignment
=====================

USDA APHIS Veterinary Services (VS) Mycobacterium tuberculosis complex, mainly M. bovis, and Brucella sp. genotyping from whole genome sequence (WGS) outputting BAM, VCF, SNP tables and phylogentic trees. 

SUMMARY
========

Using BWA alignments VCFs are created for each isolate using GATK's haplotype caller.  The overview below discribes a two step pipeline, briefly: step 1 - outputing VCFs and step 2 - gathering those VCFs to produce SNP alignments, tables and phylogenetic trees of grouped isolates.  This two step pipeline has been developed to provide routine diagnostics results rapidly, with quality assurance and easy to intrepret reporting.

## Script 1
Step 1 is fairly straight forward.  Our main workflow includes Mycobacterium tuberculosis complex and Brucella sp. and therefore the script has been optimized for such.  The script begins by selecting the best reference and determining the spoligtype for TB complex isolates and MLST for Brucella speices.  The script then goes on to using what is now a fairly standard bacterial SNP calling pipeline where BWA aligns and GATK calls SNPs.  However, in addition what is output by the HaplotypeCaller, Map zero positions are added to VCF.  This is not a standard GATK option and is done using shell commands.  Including these Map zero positions allows a more accurate SNP summary to be represented in step 2.

## Script 2
Step 2 is called on VCFs output from step 1.  References chosen in script 1 have been selected because they have been found to be relatively close to the isolate.  The closer the reference is to the isolate the less overall SNP calling error is seen.  VCFs analyzed in step 2 must all be output from the same reference.  Obviously VCFs analyzed using different references can not be used in the same comparison.


In addition to choosing a closely related reference, which minimizes SNP calling error, there are three additional external files, or dependencies, used to create high quality, informative SNP alignments.  As shown in dependency_view.jpg the reference used to build VCFs must be reflected in the three dependencies.  The three dependent files are: filter file, defining SNPs, and gbk file.

INSTALL
=======

## Python environment setup

Clone script: 


    ~$ git clone https://github.com/stuber/VCFs_to_SNP_alignment.git


Script 2 is written in Python.  Anaconda is a highly trusted Python package distrubution platform.  The script is tested using Anaconda 4.0.0.  Newer versions of Anaconda have shown to not work because some modules are incompatiable with Python >3.5.  If running a version other than 4.0.0 a new environment can be set without disrupting your current version (root) install.  See note below for installing an additional Anaconda environment.  Use Anaconda's default installation except when asked if to prepend to PATH choose yes.

Install Anaconda if not already installed.

    ~$ wget https://repo.continuum.io/archive/Anaconda3-4.0.0-Linux-x86_64.sh
    
    ~$ bash Anaconda3-4.0.0-Linux-x86_64.sh
    
Go into the python interpretter by typing $ python
Your version will be Python 3.5.1 |Anaconda 4.0.0 (64-bit).  If it is not, something did not run correctly.
Once Anaconda 4.0.0 is installed

    $ easy_install pyvcf 
    
    $ easy_install biopython
    
    $ pip install xvfbwrapper

<!--    $ easy_install ete3 -->
<!--    -->
<!--If PyQt4 is not found when trying to install ete3, install using:-->

    $ pip uninstall ete3 #uninstall

    $ conda install -c etetoolkit ete3=3.0.0b36
<!---->

<!---->
<!--    $ conda install pyqt=4-->
<!--    -->
<!--    $ easy_install ete3-->
    
Install ete3 
            
RAxML is the single program outside of the Python environment that is needed.  It must be in your PATH as: raxmlHPC-SSE3, raxmlHPC-PTHREADS-AVX2, or raxml.  In my experience installing raxmlHPC-SSE3 is the most universal but if you have the correct computer architecture running raxmlHPC-PTHREADS-AVX2 will be faster.  Below are brief instructions to install RAxML.

Download RAxML from https://github.com/stamatak/standard-RAxML.  Download to a desired location.

    ~$ git clone https://github.com/stamatak/standard-RAxML.git

    ~$ cd standard-RAxML/

    ~$ make -f Makefile.SSE3.gcc

    ~$ rm *.o

    ~$ sudo ln -s /home/user/standard-RAxML/raxmlHPC-SSE3 /usr/local/bin/raxml # create softlink in PATH

    ~$ which raxml #check that it will be found in $PATH    

## Dependency setup
By default script dependencies are expected to be in your home directory.  To install dependencies run the command below with your current working directory set to your home directory.  Check this repo periodically for updates.

~$ git clone https://github.com/stuber/dependencies.git

## Test
Use files bundled with dependencies to test.  Make working directory that containing VCFs and call script.  Call on test files with the command below.
~/dependencies/vcf_test_files $  /path/to/VCFs_to_SNP_alignment/script2.py -s bovis

###Note:  Adding an additional environment

If a new Anaconda environment is needed without making changes to your current:
        
    $ conda create -n anaconda400 anaconda=4.0.0 anaconda
    
    To activate this environment, use:
    
    > source activate anaconda400
    
    To deactivate this environment, use:
    
    > source deactivate anaconda400

OVERVIEW
========

Importance of step 1.  This is where the SNPs are called.  The final analysis is only going to be as good as the SNPs being called.  It doesn't matter which programs are used, per se, as long as each SNP called (or not called) can be justified and validated.

Because of this step 1 and step 2 are very much linked to each other.


OBJECTIVE
==========

