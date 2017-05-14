#!/usr/bin/env python

import os
import sys
import subprocess
import glob
import re
import shutil
import time
import csv
import xlrd
import vcf
import xlsxwriter
from xvfbwrapper import Xvfb
from collections import Counter
from datetime import datetime
from optparse import OptionParser
from concurrent import futures
from collections import OrderedDict
import multiprocessing
from Bio import SeqIO
import zipfile
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email import encoders

home = os.path.expanduser("~")

vdisplay = Xvfb()
vdisplay.start()

root = str(os.getcwd())
error_list = []

try:
    subprocess.call("raxml")
    sys_raxml = "raxml"
    print ("%s found" % sys_raxml)
except OSError:
    print ("looking for RAxML")
    try:
        subprocess.call("raxmlHPC-SSE3")
        sys_raxml = "raxmlHPC-SSE3"
        print ("%s found" % sys_raxml)
    except OSError:
        print ("looking for RAxML")
        try:
            subprocess.call("raxmlHPC-PTHREADS-AVX2")
            sys_raxml = "raxmlHPC-PTHREADS-AVX2"
            print ("RAxML found")
        except OSError:
            print ("#####RAxML is not in you PATH")
            print ("#####See help page for support")
            sys.exit(0)

print ("\n\nRAxML found in $PATH as: %s" % sys_raxml)

#set cpu usage
cpu_count = multiprocessing.cpu_count()

if cpu_count < 20:
    raxml_cpu = 2
else:
    raxml_cpu = int(cpu_count/10)

cpu_count = int(cpu_count/2.3)

htmlfile_name = root + "/summary_log.html"
htmlfile = open(htmlfile_name, 'at')

startTime = datetime.now()
print ("\n\n*** START ***\n")
print ("Start time: %s" % startTime)


### PARSE ARGUMENTS
parser = OptionParser()
parser.add_option('-s', '--species', action='store', dest='species', help='provide species/reference specific for VCFs, possible options: bovis, h37, ab1, ab3, suis1, suis3, suis4, mel1, mel2, mel3, canis, ceti1, ceti2', metavar='<REQUIRED options: bovis, h37, ab1, ab3, suis1, mel1, mel2, mel3, canis, ceti1, ceti2')
parser.add_option('-a', '--all_vcf', action='store_true', dest='all_vcf', help='make tree using all VCFs')
parser.add_option('-e', '--elite', action='store_true', dest='elite', help='create a tree with on elite sample representation')
parser.add_option('-d', '--debug', action='store_true', dest='debug_call', help='debug, run loops withouth future')
parser.add_option('-p', '--pilon', action='store_true', dest='pilon', help='pilon vcfs, run pilon produced vcfs')
parser.add_option('-i', '--inhouse', action='store_true', dest='inhouse', help='in-house variable setting')
parser.add_option('-u', '--upload', action='store_true', dest='upload', help='upload files to the bioinfo drive')
parser.add_option('-f', '--filter', action='store_true', dest='filter', help='Find possible positions to filter')
parser.add_option('-m', '--email', action='store', dest='email', help='email recipients: all, s, tod, jess, suelee, chris, email_address')
parser.add_option('-x', '--xserver', action='store_true', dest='xserver', help='if functional ete3 provides pdf and svg of tree')

(options, args) = parser.parse_args()
print ("SET OPTIONS: ")
print (options)

if options.email == "all":
    email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"
elif options.email == "tod":
    email_list = "tod.p.stuber@usda.gov"
elif options.email == "jess":
    email_list = "jessica.a.hicks@usda.gov"
elif options.email == "chris":
    email_list = "tod.p.stuber@usda.gov, christine.r.quance@usda.gov"
elif options.email == "suelee":
    email_list = "tod.p.stuber@usda.gov, suelee.robbe-austerman@aphis.usda.gov"
else:
    email_list = options.email

if options.species == "bovis" or options.species == "h37" or options.species == "ab1" or options.species == "ab3" or options.species == "suis1" or options.species == "suis3" or options.species == "suis4" or options.species == "mel1" or options.species == "mel2" or options.species == "mel3" or options.species == "canis" or options.species == "ceti1" or options.species == "ceti2":
    print('options.species True %s' % options.species)
else:
    print('options.species False %s' % options.species)
    print ("MUST PROVED A SPECIES")
    parser.print_help()
    sys.exit(0)

if options.all_vcf:
    print ("all_vcf %s " % options.all_vcf)
else:
    print ("all_vcf %s " % options.all_vcf)

# DIRECTORY TEST AND BACKUP
if getattr(sys, 'frozen', False):
    script_used = os.path.realpath(sys.executable)
elif __file__:
    script_used = os.path.realpath(__file__)

print ("\nScript used: %s \n" % script_used)

# make backup
os.makedirs('starting_files')
all_starting_files = glob.glob('*vcf')
for i in all_starting_files:
    shutil.copy(i, 'starting_files')

##################
# FUNCTIONS
##################

def update_directory():
    if options.inhouse:
        if os.path.isdir("/bioinfo11/TStuber/Results"):
            upload_to = "/bioinfo11/TStuber/Results"
            remote="/bioinfo11/TStuber/Results" + dependents_dir
            if os.path.isdir("/Users/Shared"):
                dep_path = "/Users/Shared"
                dir_split = dependents_dir.split('/')[1:]
                for i in dir_split:
                    dep_path += '/' + i
                    if not os.path.exists(dep_path):
                        os.makedirs(dep_path)
                local = "/Users/Shared" + dependents_dir
            elif os.path.isdir("/home/shared"):
                dep_path = "/home/shared"
                dir_split = dependents_dir.split('/')[1:]
                for i in dir_split:
                    dep_path += '/' + i
                    if not os.path.exists(dep_path):
                        os.makedirs(dep_path)
                local = "/home/shared" + dependents_dir
        elif os.path.isdir("/Volumes/root/TStuber/Results"):
            upload_to = "/Volumes/root/TStuber/Results"
            remote="/Volumes/root/TStuber/Results" + dependents_dir
            if os.path.isdir("/Users/Shared"):
                dep_path = "/Users/Shared"
                dir_split = dependents_dir.split('/')[1:]
                for i in dir_split:
                    dep_path += '/' + i
                    if not os.path.exists(dep_path):
                        os.makedirs(dep_path)
                local = "/Users/Shared" + dependents_dir
            elif os.path.isdir("/home/shared"):
                dep_path = "/home/shared"
                dir_split = dependents_dir.split('/')[1:]
                for i in dir_split:
                    dep_path += '/' + i
                    if not os.path.exists(dep_path):
                        os.makedirs(dep_path)
                local = "/home/shared" + dependents_dir
        else:
            upload_to ="not_found"
            remote = "not_found"
            if os.path.isdir("/Users/Shared" + dependents_dir):
                local = "/Users/Shared" + dependents_dir
            elif os.path.isdir("/home/shared" + dependents_dir):
                local = "/home/shared" + dependents_dir
            else:
                print ("\n\n***ERROR: CHECK PATH TO SCRIPT DEPENDENTS")
                sys.exit("SCRIPT EXITED")

        print ("\nSET LOCATIONS")
        print ("\tremote: %s " % remote)
        print ("\tlocal: %s " % local)

        if remote == "not_found":
            upload_to ="not_found"
            print ("\nUSING LOCAL SCRIPT DEPENDENTS, UNABLE TO FIND NETWORK DIRECTORY")
            return upload_to, remote, local

        else:
            print ("\nLOCALLY UPDATING SCRIPT DEPENDENTS")
            print ("\t%s --> %s" % (remote, local))
            if os.path.isdir(local):
                shutil.rmtree(local)
                shutil.copytree(remote, local)
            return upload_to, remote, local
    else:
        upload_to ="not_found"
        remote = "no remote"
        script_location = home # points to home directory
        local = script_location + "/dependencies" + dependents_dir # sets dependencies directory to home directory
        print ("Path to dependents: %s" % local)
        return upload_to, remote, local

def zip(src, dst):
    print ("\nZipping files...\n")
    zf = zipfile.ZipFile("%s.zip" % (dst), "w", zipfile.ZIP_DEFLATED)
    abs_src = os.path.abspath(src)
    for dirname, subdirs, files in os.walk(src):
        for filename in files:
            absname = os.path.abspath(os.path.join(dirname, filename))
            arcname = absname[len(abs_src) + 1:]
            zf.write(absname, arcname)
    zf.close()



# Change file names
def change_names():
    code_dictionary = {}
    try:
        wb = xlrd.open_workbook(genotypingcodes)
        ws = wb.sheet_by_index(0)
        for rownum in range(ws.nrows):
            new_name = str(ws.row_values(rownum)[0])
            new_name = new_name.rstrip()
            new_name = re.sub('[\/() ]', '_', new_name)
            new_name = re.sub('#', 'num', new_name)
            new_name = re.sub('_-', '_', new_name)
            new_name = re.sub('-_', '_', new_name)
            new_name = re.sub('__+', '_', new_name)
            new_name = re.sub('_$', '', new_name)
            new_name = re.sub('-$', '', new_name)
            new_name = re.sub(',', '', new_name)
            try:
                elite_test = ws.row_values(rownum)[1]
            except IndexError:
                #print ("except IndexError: when changing names")
                elite_test = ""
            code_dictionary.update({new_name:elite_test})
    except FileNotFoundError:
        print ("\n#### except: FileNotFoundError, there was not a \"genotypingcodes\" file given to change names\n")


def get_annotations(all_positions):
    print ("Getting annotations")
    dict_annotation = {}
    in_annotation_as_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))
    for each_absolute_pos in all_positions:
        pos_found = False
        each_absolute_pos = each_absolute_pos.split("-")
        chrom = each_absolute_pos[0]
        chrom.rstrip()
        #print ("chrom %s" % chrom)
        pos = each_absolute_pos[1]
        pos.rstrip()
        pos = int(pos)
        #print ("pos %s" % pos)
        for each_key, each_value in in_annotation_as_dict.items():
            if chrom == each_key:
                for feature in each_value.features:
                    if pos in feature and "CDS" in feature.type:
                        myproduct = "none list"
                        mylocus = "none list"
                        mygene = "none list"
                        for p in feature.qualifiers['product']:
                            myproduct = p
                        for l in feature.qualifiers['locus_tag']:
                            mylocus = l
                        if "gene" in feature.qualifiers:
                            gene = feature.qualifiers['gene']
                            for g in gene:
                                mygene = g
                        myout = myproduct + ", gene: " + mygene + ", locus_tag: " + mylocus
                        pos_found = True
        if pos_found == False:
            myout = "No annotated product"
        dict_annotation.update({chrom + "-" + str(pos):myout})
        #print ("myout %s" % myout)
    return (dict_annotation)

###################################################################
###################################################################
###################################################################

### SET PARAMETERS
if options.species == "bovis":
    
    qual_gatk_threshold = 150
    N_gatk_threshold = 200
    
    #Remove network path at and left of "Results"
    dependents_dir="/mycobacterium/tbc/tbbov/script_dependents/python/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/mycobacterium/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")
    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_002945.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/mycobacterium/tbc/tbbov/script2"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    filter_files = script_dependents + "/filter_files"
    print ("filter_files %s" % filter_files)
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:
        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if not options.email:
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

