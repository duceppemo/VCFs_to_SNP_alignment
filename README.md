VCFs to SNP alignment
=====================

USDA APHIS Veterinary Services (VS) Mycobacterium tuberculosis complex, mainly M. bovis, and Brucella sp. genotyping from whole genome sequence (WGS) outputting BAM, VCF, SNP tables and phylogentic trees. 

SUMMARY
========

Using BWA alignments VCFs are created for each isolate using GATK's haplotype caller.  The overview discribes a two step pipeline: step 1 outputing VCFs and step 2 gathering those VCFs to produce SNP alignments, tables and phylogenetic trees of grouped isolates.  This two step pipeline has been developed to provide routine diagnostics results rapidly, with quality assurance and easy to intrepret reporting.

## Script 1
See step 1 overview at:
In short, step 1 is fairly straight forward.  The main workflow includes Mycobacterium tuberculosis complex and Brucella sp. and therefore this script has been optimized for such.  It begins by selecting the best reference and determining Spoligtype for TB complex genomes and MLST for Brucella sp. genomes.  The script then goes onto using what has now a standard bacterial SNP calling pipeline where BWA aligns and GATK calls SNPs, however in addition Map zero positions are added to the HaplotypeCaller VCF, which is not a standard GATK option and is  done using shell commands.  Including these Map zero positions allows a more accurate SNP summary to be represented in step 2.

## Script 2
Step 2 is called on VCFs output from step 1.  References chosen in script 1 have been selected because they have been found to be relatively close to the isolate.  The closer the reference is to the isolate the less overall SNP calling error will be seen.  VCFs analyzed in step 2 must all be output from the same reference.  Obviously VCFs analyzed using different references can not be used in the same comparison.  Because of this step 1 and step 2 are very much linked to each other.

In addition to choosing a closely related reference, which minimizes SNP calling error, there are three additional external files, or dependencies, used to create high quality, informative SNP alignments.  As shown in dependency_view.jpg the reference used to build the VCF must be reflected in the three dependencies.

INSTALL
=======

## Python environment setup
Script 2 is written in Python.  Anaconda is a highly trusted Python package distrubution platform.  The script is tested using Anaconda 4.0.0.  Other, newer, versions have shown not to work because some modules are incompatiable with Python >3.5.  If running a version other than 4.0.0 a new environment can be set without disrupting your current version (root) install.  Use Anaconda's default installation except when asked if to prepend to PATH choose yes.

        If a new Anaconda environment is needed without making changes to your current:
        
            $ conda create -n anaconda400 anaconda=4.0.0 anaconda
            
            To activate this environment, use:
            
            > source activate anaconda400
            
            To deactivate this environment, use:
            
            > source deactivate anaconda400
            

Go into the python interpretter by typing $ python
Your version will be Python 3.5.1 |Anaconda 4.0.0 (64-bit).  If it is not, something did not run correctly.
Once Anaconda 4.0.0 is installed

    $ easy_install ete3 
    
    $ easy_install pyvcf 
    
    $ easy_install biopython
    
    $ pip install xvfbwrapper
    

## Dependency setup
By default script dependencies are expected to be in your home directory.  To install dependencies run the command below with your current working directory set to your home directory.  Check this repo periodically for updates.

~$ git clone https://github.com/stuber/dependencies.git

OVERVIEW
========

Importance of step 1.  This is where the SNPs are called.  The final analysis is only going to be as good as the SNPs being called.  It doesn't matter which programs are used, per se, as long as each SNP called (or not called) can be justified and validated.

OBJECTIVE
==========

