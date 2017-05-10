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
import pandas as pd # pandas 0.18.1 tested, 0.19.2 does not work
import numpy as np
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
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces, AttrFace

home = os.path.expanduser("~")

#os.environ["DISPLAY"]=":99"
#xvfb = subprocess.Popen(['Xvfb', ':99']) # allows not needing to use -X flag when ssh'ing into session.
#Manage headless displays with Xvfb (X virtual framebuffer
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
#        upload_to ="not_found"
#        remote = "no remote"
#        script_location = os.path.dirname(script_used) #Points to location in pyinstaller package
#        local = script_location + "/dependencies" # sets to pyinstaller dependencies directory
#        print ("Path to dependents: %s" % local)
#        return upload_to, remote, local
    ### change to dependencies residing in home directory
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

# Test for duplicate samples
def test_duplicate():
    dup_list = []
    list_of_files = glob.glob('*vcf')
    for line in list_of_files:
        line=re.sub(r'(.*)[_.].*', r'\1', line)
        dup_list.append(line)
    # find duplicates in list
    duplicates = [k for k,v in Counter(dup_list).items() if v>1]
    if len(duplicates) > 0:
        print ("Duplicates Found: %s " % duplicates)
        print ("\n***Error:  Duplicate VCFs")
        sys.exit(0)
    else:
        print ("\nno duplicate VCFs\n")

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

    names_not_changed = []
    list_of_files = glob.glob('*vcf')
    for each_vcf in list_of_files:
        vcf_found = False
        vcf_pretext = re.sub(r'(.*?)[._].*', r'\1', each_vcf) # ? was needed to make greedy, in my view the regex was searching right to left without it.
        vcf_pretext = vcf_pretext.rstrip()
        if vcf_pretext[-1].isdigit():
            myregex = re.compile(vcf_pretext + '_.*') #if number require a underscore at end (writen with 16-0338, both TB and acc number, in mind)
        else:
            myregex = re.compile(vcf_pretext + '.*') #if letter do not put a underscore at end (writen with MI in mind)
        for k, v in code_dictionary.items():
            try:
                if myregex.search(k):
                    os.rename(each_vcf, k + ".vcf")
                    vcf_found = True
            except FileNotFoundError:
                print ("except FileNotFoundError %s" % each_vcf)
        if vcf_found == False:
                    names_not_changed.append(each_vcf)
    names_not_changed = set(names_not_changed) # remove duplicates

    if options.elite:
        list_of_files = []
        list_of_files = glob.glob('*vcf')
        if not os.path.exists("temp_hold"):
            print ("making temp_hold directory")
            os.makedirs("temp_hold") # make all_vcfs if none exists
        for each_vcf in list_of_files:
            time_test = time.time() - os.path.getmtime(each_vcf) < (1 * 24 * 60 *60) # 1day * (24*60*60)sec in day
            print ("%s each_vcf" % each_vcf)
            vcf_pretext = re.sub(r'(.*?)[._].*', r'\1', each_vcf) # ? was needed to make greedy, in my view the regex was searching right to left without it.
            vcf_pretext = vcf_pretext.rstrip()
            myregex = re.compile(vcf_pretext + '.*')
            if time_test:
                print ("time_test true %s" % each_vcf)
                shutil.copy(each_vcf, "temp_hold")
            else:
                for k, v in code_dictionary.items():
                    if myregex.search(k):
                        try:
                            print ("##### %s" % time_test)
                            if v == "Yes": # if marked yes in column 2 of genotyping codes
                                print ("marked yes %s" % each_vcf)
                                shutil.copy(each_vcf, "temp_hold") # if "Yes" then moved to temp_hold
                            else:
                                print ("file will be discarded %s" % each_vcf)
                        except FileNotFoundError:
                            print ("except FileNotFoundError %s" % each_vcf)
                os.remove(each_vcf)
        shutil.rmtree('starting_files')
        os.makedirs('starting_files')
        os.renames('temp_hold', 'starting_files')
        list_of_files = glob.glob('starting_files/*vcf')
        file_number = len(list_of_files) # update the file_number to present on summary
        for each_vcf in list_of_files:
            shutil.copy(each_vcf, root)
        all_starting_files = glob.glob('*vcf')
        print (file_number)

    return names_not_changed

# Get filters set up
def get_filters(excelinfile, filter_files):
    for i in glob.glob(filter_files + "/*"):
        os.remove(i)

    wb = xlrd.open_workbook(excelinfile)
    sheets = wb.sheet_names()
    for sheet in sheets:
        ws = wb.sheet_by_name(sheet)

        myrange = lambda start, end: range(start, end+1)

        for colnum in range(ws.ncols): # for each column in worksheet
            file_out = filter_files + "/" + ws.col_values(colnum)[0] + ".txt" # column header naming file
            write_out = open (file_out, 'at')
            mylist = ws.col_values(colnum)[1:] # list of each field in column, minus the header
            mylist = [x for x in mylist if x] # remove blank cells
            for value in mylist:
                value = str(value)
                value = value.replace(sheet + "-", '')
                if "-" not in value:
                    value=int(float(value)) # change str to float to int
                    print (sheet + "-" + str(value), file=write_out)
                elif "-" in value:
                    value = value.split("-")
                    for i in range(int(value[0]), int(value[1]) + 1 ):
                        print (sheet + "-" + str(i), file=write_out)
    write_out.close()

# Group files
def group_files(each_vcf):
    list_pass = []
    list_amb = []
    dict_amb = {}
    malformed = []
    group_calls = []
    passing = True
    
    ###
    # Fix common VCF errors
    if options.debug_call:
        print ("FIXING FILE: " + each_vcf)
    temp_file = each_vcf + ".temp"
    write_out=open(temp_file, 'w') #r+ used for reading and writing to the same file

    with open(each_vcf, 'r') as file:
        try:
            for line in file:
                if line.rstrip(): # true if not empty line'^$'
                    line = line.rstrip() #remove right white space
                    line = re.sub('"AC=', 'AC=', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('"$', '', line)
                    line = re.sub('GQ:PL\t"', 'GQ:PL\t', line)
                    line = re.sub('[0-9]+\tGT\t.\/.$', '999\tGT:AD:DP:GQ:PL\t1/1:0,80:80:99:2352,239,0', line)
                    line = re.sub('^"', '', line)
                    if line.startswith('##'):
                        line = line.split('\t')
                        line = ''.join(line[0])
                    if not line.startswith('##'):
                        line = re.sub('"', '', line)
                        line = line.split('\t')
                        line = "\t".join(line[0:10])
                        print(line, file=write_out)
                    else:
                        print(line, file=write_out)
        except IndexError:
            print ("##### IndexError: Deleting corrupt VCF file: " + each_vcf)
            malformed = "##### IndexError: Deleting corrupt VCF file: " + each_vcf
            os.remove(each_vcf)
        except UnicodeDecodeError:
            print ("##### UnicodeDecodeError: Deleting corrupt VCF file: " + each_vcf)
            malformed = "##### UnicodeDecodeError: Deleting corrupt VCF file: " + each_vcf
            os.remove(each_vcf)
            
    write_out.close()
    os.rename(temp_file, each_vcf)
    ###
    
    try:
        vcf_reader = vcf.Reader(open(each_vcf, 'r'))
        ### PUT VCF NAME INTO LIST, capturing for htmlfile
        group_calls.append(each_vcf)
            # for each single vcf getting passing position
        for record in vcf_reader:
            chrom = record.CHROM
            position = record.POS
            absolute_positon = str(chrom) + "-" + str(position)
            # find quality SNPs and put absolute positions into list
            if options.pilon:
                filter=record.FILTER
                # find quality SNPs and put absolute positions into list
                try:
                    record_alt_length = len(record.ALT[0])
                except TypeError:
                    record_alt_length = 0
                try:
                    record_ref_length = len(record.REF)
                except TypeError:
                    record_alt_length = 0
                
                if not filter and record_alt_length == 1 and record_ref_length == 1 and record.QUAL > 50:
                    list_pass.append(absolute_positon)
                # capture ambigous defining SNPs in htmlfile
                elif record.FILTER == ['Amb']:
                    list_amb.append(absolute_positon)
            else:
                try:
                    record_alt_length = len(record.ALT[0])
                except TypeError:
                    record_alt_length = 0
                try:
                    record_ref_length = len(record.REF)
                except TypeError:
                    record_alt_length = 0
                try:
                    if str(record.ALT[0]) != "None" and record_ref_length == 1 and record_alt_length == 1 and record.INFO['AC'][0] == 2 and record.QUAL > qual_gatk_threshold and record.INFO['MQ'] > 45:
                        list_pass.append(absolute_positon)
                    # capture ambigous defining SNPs in htmlfile
                    elif str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 1:
                        list_amb.append(absolute_positon)
                except ZeroDivisionError:
                    print ("bad line in %s at %s" % (each_vcf, absolute_positon))

        for key in inverted_position.keys():
            if key not in list_pass:
                print ("key %s not in list_pass" % key)
                directory = inverted_position[key]
                print("*** INVERTED POSITION FOUND *** PASSING POSTION FOUND: \t%s\t\t%s" % (each_vcf, directory))
                if not os.path.exists(directory):
                    try:
                        os.makedirs(directory)
                    except FileExistsError:
                        null = "null"
                shutil.copy(each_vcf, directory)
                ### ADD GROUP TO LIST
                group_calls.append(directory)

        #if passing:
        # if a passing position is in the defining SNPs
        for passing_position in list_pass:
            # normal grouping
            if passing_position in defining_snps:
                directory = defining_snps[passing_position]
                print("PASSING POSTION FOUND: \t%s\t\t%s" % (each_vcf, directory))
                if not os.path.exists(directory):
                    try:
                        os.makedirs(directory)
                    except FileExistsError:
                        null = "null"
                shutil.copy(each_vcf, directory)
                ### ADD GROUP TO LIST
                group_calls.append(directory)
                
        # find mixed isolates if defining snp is ambigous
        for amb_position in list_amb:
            if amb_position in defining_snps:
                directory = defining_snps[amb_position]
                dict_amb.update({each_vcf + "\t" + directory:amb_position})
                ### ADD AMBIGIOUS CALL TO LIST
                group_calls.append("*" + directory + "-mix")
        # if -a or -e (non elites already deleted from the analysis) copy all vcfs to All_VCFs
        if options.all_vcf or options.elite:
            if not os.path.exists("All_VCFs"):
                os.makedirs("All_VCFs")
            shutil.move(each_vcf, "All_VCFs")
        else:
            try:
                os.remove(each_vcf)
            except FileNotFoundError:
                print ("file deleted: %s " % each_vcf)
                malformed.append(each_vcf)
        #print (dict_amb, group_calls, malformed)
        
        try:
            some_object_iterator = iter(group_calls)
        except TypeError:
            group_calls = []

        try:
            some_object_iterator = iter(dict_amb)
        except TypeError:
            dict_amb = {}

        try:
            some_object_iterator = iter(malformed)
        except TypeError:
            os.remove(each_vcf)
            print ("TypeError: corrupt VCF, removed %s " % each_vcf)
            malformed = "TypeError: corrupt VCF, removed %s " % each_vcf

    except ZeroDivisionError:
        os.remove(each_vcf)
        print ("ZeroDivisionError: corrupt VCF, removed %s " % each_vcf)
        malformed = "ZeroDivisionError: corrupt VCF, removed %s " % each_vcf
    except ValueError:
        os.remove(each_vcf)
        print ("ValueError: corrupt VCF, removed %s " % each_vcf)
        malformed = "ValueError: corrupt VCF, removed %s " % each_vcf
    except UnboundLocalError:
        os.remove(each_vcf)
        print ("UnboundLocalError: corrupt VCF, removed %s " % each_vcf)
        malformed = "UnboundLocalError: corrupt VCF, removed %s " % each_vcf
    except TypeError:
        os.remove(each_vcf)
        print ("TypeError: corrupt VCF, removed %s " % each_vcf)
        malformed = "TypeError: corrupt VCF, removed %s " % each_vcf
    except SyntaxError:
        os.remove(each_vcf)
        print ("SyntaxError: corrupt VCF, removed %s " % each_vcf)
        malformed = "SyntaxError: corrupt VCF, removed %s " % each_vcf
    except KeyError:
        os.remove(each_vcf)
        print ("KeyError: corrupt VCF, removed %s " % each_vcf)
        malformed = "KeyError: corrupt VCF, removed %s " % each_vcf
    except StopIteration:
        print ("StopIteration: %s" % each_vcf)
        malformed = "KeyError: corrupt VCF, removed %s " % each_vcf

    a = group_calls[0:1]
    b = sorted(group_calls[1:]) # order the groups
    for i in b:
        a.append(i) # a is group_calls
        group_calls = a
    return dict_amb, group_calls, malformed

# Table to Excel file
def excelwriter(filename):
    orginal_name=filename
    filename = filename.replace(".txt",".xlsx")
    wb = xlsxwriter.Workbook(filename)
    ws = wb.add_worksheet("Sheet1")
    with open(orginal_name,'r') as csvfile:
        table = csv.reader(csvfile, delimiter='\t')
        i = 0
        for row in table:
            ws.write_row(i, 0, row)
            i += 1

    col = len(row)
    col = col + 1
    #print (i, "x", col)

    formatA = wb.add_format({'bg_color':'#58FA82'})
    formatG = wb.add_format({'bg_color':'#F7FE2E'})
    formatC = wb.add_format({'bg_color':'#0000FF'})
    formatT = wb.add_format({'bg_color':'#FF0000'})
    formatnormal = wb.add_format({'bg_color':'#FDFEFE'})
    formatlowqual = wb.add_format({'font_color':'#C70039', 'bg_color':'#E2CFDD'})
    formathighqual = wb.add_format({'font_color':'#000000', 'bg_color':'#FDFEFE'})
    formatambigous = wb.add_format({'font_color':'#C70039', 'bg_color':'#E2CFDD'})
    formatN = wb.add_format({'bg_color':'#E2CFDD'})

    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                          'criteria':'containing',
                          'value':60,
                          'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':59,
                        'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':58,
                        'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':57,
                        'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':56,
                        'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':55,
                        'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':54,
                        'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':53,
                        'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':52,
                        'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':51,
                        'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':50,
                        'format':formathighqual})
    ws.conditional_format(i-2,1,i-2,col-1, {'type':'text',
                        'criteria':'not containing',
                        'value':100,
                        'format':formatlowqual})

    ws.conditional_format(2,1,i-3,col-1, {'type':'cell',
                        'criteria':'==',
                        'value':'B$2',
                        'format':formatnormal})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'A',
                        'format':formatA})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'G',
                        'format':formatG})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'C',
                        'format':formatC})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'T',
                        'format':formatT})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'S',
                        'format':formatambigous})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'Y',
                        'format':formatambigous})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'R',
                        'format':formatambigous})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'W',
                        'format':formatambigous})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'K',
                        'format':formatambigous})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'M',
                        'format':formatambigous})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'N',
                        'format':formatN})
    ws.conditional_format(2,1,i-3,col-1, {'type':'text',
                        'criteria':'containing',
                        'value':'-',
                        'format':formatN})

    ws.set_column(0, 0, 30)
    ws.set_column(1, col-1, 2)
    ws.freeze_panes(2, 1)
    format_rotation = wb.add_format({'rotation':'90'})
    ws.set_row(0, 140, format_rotation)
    formatannotation = wb.add_format({'font_color':'#0A028C', 'rotation':'-90', 'align':'top'})
    #set last row
    ws.set_row(i-1, 400, formatannotation)

    wb.close()

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

def find_positions(filename):
    found_positions = {}
    vcf_reader = vcf.Reader(open(filename, 'r'))
    try:
        for record in vcf_reader:
            chrom = record.CHROM
            position = record.POS
            absolute_positon = str(chrom) + "-" + str(position)
            filter=record.FILTER
            
            # Usable positins are those that:

            # ADD PARAMETERS HERE TO CHANGE WHAT'S SNP WILL BE USED
            # IF NOT FOUND HERE THE SNP WILL BE IGNORED.  WILL NOT BE REPRESENTED.  HARD REMOVAL
            
            ## GATK parameters
            # str(record.ALT[0]) != "None" --> filter deletions
            # len(record.REF) == 1 --> filter bad ref call with 2 nt present
            # len(record.ALT[0]) == 1 --> filter bad alt call with 2 nt present
            # record.heterozygosity == 0.0 --> filter AC=1, heterozygosity.
            # record.QUAL > 150 --> filter poor quality
            # record.INFO['MQ'] --> filter low map quality
            if options.pilon:
                if not filter and len(record.REF) == 1 and len(record.ALT[0]) == 1 and record.QUAL > 50 and record.INFO['MQ'] > 30:
                    found_positions.update({absolute_positon:record.REF})
            else:
                try:
                    if str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 2 and len(record.REF) == 1 and len(record.REF) == 1 and record.QUAL > qual_gatk_threshold:
                        found_positions.update({absolute_positon:record.REF})
                except KeyError:
                    pass
    except ZeroDivisionError:
        print ("ZeroDivisionError error found")
    except ValueError:
        print ("ValueError error found")
    except UnboundLocalError:
        print ("UnboundLocalError error found")
    except TypeError:
        print ("TypeError error found")

    return found_positions

def sort_table(intable, ordered, out_org):
    mytable = pd.read_csv(intable, sep='\t')
    mytable=mytable.set_index('reference_pos')

    # order list is from tree file
    # gives order for samples to be listed in table to be phylogenetically correct
    ordered_list = []
    with open(ordered) as infile:
        for i in infile:
            i = i.rstrip()
            ordered_list.append(i)
    # sinces this is set as the mytable index do not include in ordering
    ordered_list.remove('reference_pos')

    # reorder table based on order of list
    mytable = mytable.reindex(ordered_list)

    # count number of SNPs in each column
    snp_per_column = []
    for column_header in mytable:
        count = 0
        column = mytable[column_header]
        # for each element in the column
        for element in column:
            if element != column[0]:
                count = count + 1
        snp_per_column.append(count)
        #print ("the count is: %s" % count)
    row1 = pd.Series (snp_per_column, mytable.columns, name="snp_per_column")

    # get the snp count per column
    # for each column in the table
    snp_from_top = []
    for column_header in mytable:
        count = 0
        column = mytable[column_header]
        # for each element in the column
        # skip the first element
        for element in column[1:]:
            if element == column[0]:
                count = count + 1
            else:
                break
        snp_from_top.append(count)
        #print ("the count is: %s" % count)
    row2 = pd.Series (snp_from_top, mytable.columns, name="snp_from_top")

    mytable = mytable.append([row1])
    mytable = mytable.append([row2])

    mytable = mytable.T
    mytable = mytable.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
    mytable = mytable.T

    # remove snp_per_column and snp_from_top rows
    mytable = mytable[:-2]
    mytable.to_csv(out_org, sep='\t')

def find_filter_dict(each_vcf):
    dict_qual = {}
    dict_map = {}
    vcf_reader = vcf.Reader(open(each_vcf, 'r'))
    for record in vcf_reader:
        absolute_positon = str(record.CHROM) + "-" + str(record.POS)
        if record.QUAL:
            returned_qual = []
            returned_qual.append(record.QUAL)
        try:
            returned_map = []
            returned_map.append(record.INFO['MQ'])
        except KeyError:
            pass
        
        dict_qual[absolute_positon] = returned_qual
        dict_map[absolute_positon] = returned_qual
    return dict_qual, dict_map

def get_snps(directory):
    os.chdir(root + "/" + directory)
    print ("\n----------------------------")
    print ("\nworking on: %s " % directory)
    outdir=str(os.getcwd()) + "/"
    # FILTER position all list
    list_filter_files = glob.glob(filter_files + '/*')

    filter_file = "empty" # if filter an all_vcf file not found mark as empty
    filter_group = "empty" # if a group specific filter file is not found mark as empty
    for i in list_filter_files:
        if "-All.txt" in i:
            filter_file = i
    
    for i in list_filter_files:
        if directory  + ".txt" in i:
            filter_group = i
    
    print ("%s --> filter_file %s " % (directory, filter_file))
    print ("%s --> filter_group %s " % (directory, filter_group))
    print ("%s --> outdir %s " % (directory, outdir))
    
    files = glob.glob('*vcf')
    all_positions = {}
    
    if options.debug_call:
        for i in files:
            found_positions = find_positions(i)
            all_positions.update(found_positions)
    else:
        with futures.ProcessPoolExecutor(max_workers=cpu_count) as pool:
            for found_positions in pool.map(find_positions, files):
                all_positions.update(found_positions)

    print ("Directory %s found positions %s" % (directory, len(all_positions)))
    presize=len(all_positions)

    # Filter applied to all positions
    if not filter_file is "empty":
        with open(filter_file, 'rt') as f:
            filter_list = f.read().splitlines() #removes \n
        for pos in filter_list:
            all_positions.pop(pos, None)
        f.close()

    # Filter applied to group
    if not filter_group is "empty":
        with open(filter_group, 'rt') as f:
            filter_list = f.read().splitlines() #removes \n
        for pos in filter_list:
            all_positions.pop(pos, None)
        f.close()
    
    print ("\nDirectory: ", directory)
    print ("Total positions found: %s" % format(presize, ",d"))
    print ("Possible positions filtered %s" % format(len(filter_list), ",d"))
    print ("Positions after filtering %s\n" % format(len(all_positions), ",d"))

#####
    # NEEDS TO BE FIXED
    if options.filter:
        #write to files
        positions_to_filter = "positions_to_filter.txt"
        positions_to_filter_details = "positions_to_filter_details.txt"
        write_out_positions=open(positions_to_filter, 'w')
        write_out_details=open(positions_to_filter_details, 'w')

        files = glob.glob('*vcf')

        #calculate mean/max qual and map at all possible positions
        from collections import defaultdict
        dd_qual = {}
        dd_map = {}
        if options.debug_call:
            for each_vcf in files:
                print ("working on: %s" % each_vcf)
                dict_qual, dict_map = find_filter_dict(each_vcf)
                keys = set(dd_qual).union(dict_qual)
                no = []
                dd_qual = dict((k, dd_qual.get(k, no) + dict_qual.get(k, no)) for k in keys)
                keys = set(dd_map).union(dict_map)
                no = []
                dd_map = dict((k, dd_map.get(k, no) + dict_map.get(k, no)) for k in keys)
        else:
            with futures.ProcessPoolExecutor() as pool:
                for dict_qual, dict_map in pool.map(find_filter_dict, files):
                    keys = set(dd_qual).union(dict_qual)
                    no = []
                    dd_qual = dict((k, dd_qual.get(k, no) + dict_qual.get(k, no)) for k in keys)
                    keys = set(dd_map).union(dict_map)
                    no = []
                    dd_map = dict((k, dd_map.get(k, no) + dict_map.get(k, no)) for k in keys)

        #dict_qual=dict((k, v) for k, v in dict_qual.items() if v)
        #dict_map=dict((k, v) for k, v in dict_map.items() if v)

        ave_qual = {}
        max_qual = {}
        for k, v in dd_qual.items():
            #only use if > 3 positions have been called
            if len(v) > 3:
                ave_qual[k]=np.mean(v)
                max_qual[k]=np.max(v)

        ave_map = {}
        max_map = {}
        for k, v in dd_map.items():
            if len(v) > 3:
                ave_map[k]=np.mean(v)
                max_map[k]=np.max(v)		

        # get all possible used positions
        all_maybe_filter = []
        for k in ave_qual.keys():
            all_maybe_filter.append(k)
        for k in max_qual.keys():
            all_maybe_filter.append(k)
        for k in ave_map.keys():
            all_maybe_filter.append(k)
        for k in max_map.keys():
            all_maybe_filter.append(k)
            # remove duplicates
            all_maybe_filter = list(set(all_maybe_filter))

        #remove those in filter list
        #Filter applied to all positions
        if not filter_file is "empty":
            with open(filter_file, 'rt') as f:
                filter_list = f.read().splitlines() #removes \n
                try:
                    for pos in filter_list:
                        all_maybe_filter.pop(pos)
                except TypeError:
                    pass
                except KeyError:
                    pass
            f.close()

        # Filter applied to group
        if not filter_group is "empty":
            with open(filter_group, 'rt') as f:
                filter_list = f.read().splitlines() #removes \n
                try:
                    for pos in filter_list:
                        all_maybe_filter.pop(pos)
                except TypeError:
                    pass
                except KeyError:
                    pass
            f.close()

        # for each possible posible position check if to filter.
        for absolute_positon in all_maybe_filter:
            ave_qual_value = ave_qual[absolute_positon]
            max_qual_value = max_qual[absolute_positon]
            ave_map_value = ave_map[absolute_positon]
            max_map_value = max_map[absolute_positon]
            print ("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value))
            if max_qual_value < 1300 and ave_qual_value < 800 and max_map_value < 58 and ave_map_value < 58:
                print ("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value), file=write_out_details)
                print (absolute_positon, file=write_out_positions)
        write_out_positions.close()
        write_out_details.close()

######
    if mygbk:
        dict_annotation = get_annotations(all_positions)
        write_out=open('annotations.txt', 'w+')
        print ('reference_pos\tannotations', file=write_out)
        for k, v in dict_annotation.items():
            print ('%s\t%s' % (k, v), file=write_out)
        write_out.close()
    #sys.exit(0)

    ########
    #^
    #|
    #possible SNPs found, all positions in dictionary

    # NOW COLLECTING FOR EACH INDIVIDUAL SAMPLE/VCF
    #select positions for each sample

    #?? how to write table directly to pd dataframe ??

    out_table= outdir + directory + "-table.txt"
    table=open(out_table, 'wt')

    # write absolute positions to table
    # order before adding to file to match with ordering of individual samples below
    all_positions=OrderedDict(sorted(all_positions.items()))
    print ("reference_pos", end="\t", file=table)
    for k, v in all_positions.items():
        print(k, end="\t", file=table)
    print ("", file=table)

    list_of_files = glob.glob('*vcf')
    
    # for each vcf
    all_map_qualities={}
    for file_name in list_of_files:
        sample_map_qualities={}
        just_name = file_name.replace('.vcf', '')
        just_name = re.sub('\..*', '*', just_name) # if after the .vcf is removed there is stilll a "." in the name it is assumed the name did not get changed
        print(just_name, end="\t", file=table)
        # for each line in vcf
        vcf_reader = vcf.Reader(open(file_name, 'r'))
        sample_dict = {}
        for record in vcf_reader:
            record_position = str(record.CHROM) + "-" + str(record.POS)
            if record_position in all_positions:
                #print ("############, %s, %s" % (file_name, record_position))
                # NOT SURE THIS IS THE BEST PLACE TO CAPTURE MQ AVERAGE
                # MAY BE FASTER AFTER PARSIMONY SNPS ARE DECIDED, BUT THEN IT WILL REQUIRE OPENING THE FILES AGAIN.
                if str(record.ALT[0]) != "None" and str(record.INFO['MQ']) != "nan": #on rare occassions MQ gets called "NaN" thus passing a string when a number is expected when calculating average.
                    #print ("getting map quality:    %s          %s      %s" % (record.INFO['MQ'], file_name, str(record.POS)))
                    sample_map_qualities.update({record_position:record.INFO['MQ']})
                # ADD PARAMETERS HERE TO CHANGE WHAT'S EACH VCF REPRESENTS.
                # SNP IS REPRESENTED IN TABLE, NOW HOW WILL THE VCF REPRESENT THE CALLED POSITION
                # str(record.ALT[0]) != "None", which means a deletion as ALT
                # not record.FILTER, or rather PASSED.
                
                # check record.QUAL
                # In GATK VCFs "!= None" not used.
                if str(record.ALT[0]) != "None" and len(record.ALT[0]) == 1 and record.INFO['AC'][0] == 2:
                    sample_dict.update({record_position:record.ALT[0]})
                # same as above but take into account Ambiguious call
                #elif str(record.ALT[0]) != "None" and len(record.ALT[0]) == 1 and record.INFO['AC'][0] == 1 and record.QUAL >= N_gatk_threshold:
                elif str(record.ALT[0]) != "None" and len(record.ALT[0]) == 1 and record.INFO['AC'][0] == 1:
                    ref_alt = str(record.ALT[0]) + str(record.REF[0])
                    if ref_alt == "AG":
                        sample_dict.update({record_position:"R"})
                    elif ref_alt == "CT":
                        sample_dict.update({record_position:"Y"})
                    elif ref_alt == "GC":
                        sample_dict.update({record_position:"S"})
                    elif ref_alt == "AT":
                        sample_dict.update({record_position:"W"})
                    elif ref_alt == "GT":
                        sample_dict.update({record_position:"K"})
                    elif ref_alt == "AC":
                        sample_dict.update({record_position:"M"})
                    elif ref_alt == "GA":
                        sample_dict.update({record_position:"R"})
                    elif ref_alt == "TC":
                        sample_dict.update({record_position:"Y"})
                    elif ref_alt == "CG":
                        sample_dict.update({record_position:"S"})
                    elif ref_alt == "TA":
                        sample_dict.update({record_position:"W"})
                    elif ref_alt == "TG":
                        sample_dict.update({record_position:"K"})
                    elif ref_alt == "CA":
                        sample_dict.update({record_position:"M"})
                    else:
                        sample_dict.update({record_position:"N"})
                    # Poor calls
                elif str(record.ALT[0]) != "None" and record.QUAL < N_gatk_threshold:
                    sample_dict.update({record_position:"N"})
                # same as above but take into account Deletion call
                elif str(record.ALT[0]) == "None":
                    sample_dict.update({record_position:"-"})

        # After iterating through VCF combine dict to nested dict
        all_map_qualities.update({just_name: sample_map_qualities})

        # merge dictionaries and order
        merge_dict={}
        merge_dict.update(all_positions)
        merge_dict.update(sample_dict)
        merge_dict=OrderedDict(sorted(merge_dict.items()))
        for k, v in merge_dict.items():
            #print ("k %s, v %s" % (k, v))
            print (str(v) + "\t", file=table, end="")
        print ("", file=table)
    table.close()
    
    ## Select parsimony informative SNPs
    mytable = pd.read_csv(out_table, sep='\t')
    mytable = mytable.set_index('reference_pos')

    # drop NaN rows and columns
    mytable=mytable.dropna(axis=1)

    # SELECT PARISOMONY INFORMATIVE SNPSs
    # removes columns where all fields are the same
    parsimony=mytable.loc[:, (mytable != mytable.ix[0]).any()]
    parsimony_positions=list(parsimony)

    write_out=open("each_vcf-poslist.txt", 'wt')
    for i in parsimony_positions:
        write_out.write(i + "\n")
    write_out.close()

    parsimony.to_csv(out_table, sep="\t", index_label='reference_pos')
    table=open(out_table, 'a')

    # added corresponding reference to parsimony table
    print ("reference_call", end="\t", file=table)
    all_positions_list=list(all_positions)
    for l in parsimony_positions:
        print(all_positions.get(l), end="\t", file=table)
    print ("", file=table)
    table.close()

    # fix end, NEED TO REWRITE, WRITING TO NEW FILE IS CLUNKY FIX
    write_table= outdir + directory + "-write.txt"
    write_out=open(write_table, 'wt')
    with open(out_table, 'rt') as f:
        for line in f:
            line = line.replace('\t\n', '\n')
            print (line, file=write_out)
    write_out.close()

    #Print out fasta alignment file from table
    alignment_file= outdir + directory + ".fasta"
    write_out=open(alignment_file, 'wt')
    with open(out_table, 'rt') as f:
        count=0
        for line in f:
            if count > 0:
                line=re.sub('^', '>', line)
                line=line.replace('reference_call', 'root')
                line=line.replace('\t', '\n', 1)
                line=line.replace('\t', '')
                print (line, end="", file=write_out)
            count = count + 1
    write_out.close()

    mytable = pd.read_csv(write_table, sep='\t')
    mytable = mytable.set_index('reference_pos')
    mytable=mytable.dropna(axis=1)

    os.remove(write_table)

    # move reference to top row
    myref=mytable[-1:]
    myother=mytable[:-1]
    frames = [myref, myother]
    mytable=pd.concat(frames)
    mytable.to_csv(out_table, sep="\t", index_label='reference_pos')

    print ("\n%s table dimensions: %s" % (directory, str(mytable.shape)))

    print ("%s RAxML running..." % directory)
    rooted_tree = outdir + directory + "-rooted.tre"
    os.system("{} -s {} -n raxml -m GTRCATI -o root -p 12345 -T {} > /dev/null 2>&1" .format(sys_raxml, alignment_file, raxml_cpu))

    try:
        ordered_list_from_tree = outdir + directory + "-cleanedAlignment.txt"
        write_out=open(ordered_list_from_tree, 'w+')
        print ("reference_pos", file=write_out)
        print ("reference_call", file=write_out)
        if os.path.isfile("RAxML_bestTree.raxml"):
            with open("RAxML_bestTree.raxml", 'rt') as f:
                for line in f:
                    line=re.sub('[:,]', '\n', line)
                    line=re.sub('[)(]', '', line)
                    line=re.sub('[0-9]\.[0-9].*\n', '', line)
                    line=re.sub('root\n', '', line)
                    write_out.write(line)
            os.rename("RAxML_bestTree.raxml", "RAxML_bestTree.raxml.tre")
            write_out.close()

        out_org = outdir + directory + "-organized-table.txt"

        sort_table(out_table, ordered_list_from_tree, out_org) #function

        print ("%s Getting map quality..." % directory)
        average=lambda x: x.mean()
        all_map_qualities=pd.DataFrame(all_map_qualities)
        #ave_mq = Type: Series
        ave_mq = all_map_qualities.apply(average, axis=1)
        ave_mq = ave_mq.astype(int)
        ave_mq.to_csv('outfile.txt', sep='\t') # write to csv

        write_out=open('map_quality.txt', 'w+')
        print ('reference_pos\tmap-quality', file=write_out)
        with open('outfile.txt', 'rt') as f:
            for line in f:
                write_out.write(line)
        write_out.close()
        #os.remove('outfile.txt')

        #add_map_qualities() #***FUNCTION CALL
        #seemed pooling did not like a function with no parameters given
        quality = pd.read_csv('map_quality.txt', sep='\t')

        mytable = pd.read_csv(out_table, sep='\t')
        mytable=mytable.set_index('reference_pos')

        # order list is from tree file
        # gives order for samples to be listed in table to be phylogenetically correct
        ordered_list = []
        with open(ordered_list_from_tree) as infile:
            for i in infile:
                i = i.rstrip()
                ordered_list.append(i)
        # sinces this is set as the mytable index do not include in ordering
        ordered_list.remove('reference_pos')

        # reorder table based on order of list
        mytable = mytable.reindex(ordered_list)
        mytable.to_csv(out_table, sep='\t')

        out_sort=str(os.getcwd()) + "/" + directory + "-sorted-table.txt" #sorted
        mytable_sort = pd.read_csv(out_table, sep='\t') #sorted
        mytable_sort = mytable_sort.set_index('reference_pos') #sorted
        mytable_sort = mytable_sort.transpose() #sort
        mytable_sort.to_csv(out_sort, sep='\t', index_label='reference_pos') #sort

        out_org=str(os.getcwd()) + "/" + directory + "-organized-table.txt" #org
        mytable = pd.read_csv(out_org, sep='\t') #org
        mytable = mytable.set_index('reference_pos') #org
        mytable = mytable.transpose() #org
        mytable.to_csv(out_org, sep='\t', index_label='reference_pos') #org

        if mygbk:
            print ("%s gbk is present, getting annotation...\n" % directory)
            annotations = pd.read_csv('annotations.txt', sep='\t') #sort
            mytable_sort = pd.read_csv(out_sort, sep='\t') #sort
            mytable_sort = mytable_sort.merge(quality, on='reference_pos', how='inner') #sort
            mytable_sort = mytable_sort.merge(annotations, on='reference_pos', how='inner') #sort
            mytable_sort = mytable_sort.set_index('reference_pos') #sort
            mytable_sort = mytable_sort.transpose() #sort
            mytable_sort.to_csv(out_sort, sep='\t', index_label='reference_pos') #sort
    
            #annotations = pd.read_csv('annotations.txt', sep='\t') #org
            mytable = pd.read_csv(out_org, sep='\t') #org
            mytable = mytable.merge(quality, on='reference_pos', how='inner') #org
            mytable = mytable.merge(annotations, on='reference_pos', how='inner') #org
            mytable = mytable.set_index('reference_pos') #org
            mytable = mytable.transpose() #org
            mytable.to_csv(out_org, sep='\t', index_label='reference_pos') #org

        else:
            print ("No gbk file or no table to annotate")
            mytable_sort = pd.read_csv(out_sort, sep='\t') #sort
            mytable_sort = mytable_sort.merge(quality, on='reference_pos', how='inner') #sort
            mytable_sort = mytable_sort.set_index('reference_pos') #sort
            mytable_sort = mytable_sort.transpose() #sort
            mytable_sort.to_csv(out_sort, sep='\t', index_label='reference_pos') #sort
            # add when no annotation
            with open(out_sort, 'rt') as f:
                line=f.readline()
            f.close()
            column_count = line.count('\t') #sort
            column_count = column_count - 1 #sort
            #print ("column_count: %s" % column_count)
            with open(out_sort, 'at') as f:
                print ("no_annotation", end = '', file=f)
                print ('\t' * column_count, file=f)
            f.close()

            print ("No gbk file or no table to annotate")
            mytable = pd.read_csv(out_org, sep='\t') #org
            mytable = mytable.merge(quality, on='reference_pos', how='inner') #org
            mytable = mytable.set_index('reference_pos') #org
            mytable = mytable.transpose() #org
            mytable.to_csv(out_org, sep='\t', index_label='reference_pos') #org
            # add when no annotation
            with open(out_org, 'rt') as f:
                line=f.readline()
            f.close()
            column_count = line.count('\t')
            column_count = column_count - 1
            #print ("column_count: %s" % column_count)
            with open(out_org, 'at') as f:
                print ("no_annotation", end = '', file=f)
                print ('\t' * column_count, file=f)
            f.close()
        
        excelwriter(out_sort) #***FUNCTION CALL #sort
        excelwriter(out_org) #***FUNCTION CALL #org

        for r in glob.glob('*vcf'):
            os.remove(r)
    
    except ValueError:
        print ("##### ValueError: %s #####" % file_name)
        return
    try:
        os.remove(ordered_list_from_tree)
        os.remove('each_vcf-poslist.txt')
        os.remove('map_quality.txt')
        if mygbk:
            os.remove("annotations.txt")
        os.remove("outfile.txt")
        os.remove(out_sort)
        os.remove(out_org) # organized.txt table
        os.remove(out_table) # unorganized table
        os.remove('RAxML_info.raxml')
        os.remove('RAxML_log.raxml')
        os.remove('RAxML_parsimonyTree.raxml')
        os.remove('RAxML_result.raxml')
        os.remove(directory + '.fasta.reduced')

    except FileNotFoundError:
        pass

    ### PANDA NOTES ###
    # get the index: mytable.index
    # get columns: mytable.columns
    # get a column: mytable.AF2122_NC002945_105651, shows index (sample names)
    # get a row: mytable.ix['reference'], shows columns (positions and SNPs)
    # values: mytable.values, SNPs - series
    # strip off the bottom row: mytable[:-1]
    # get the bottom row: mytable[-1:]

    # ete3 used to make svg and pdf from trees
    # Anaconda 4.0 is needed to install ete3.  Shown to work with Anaconda 4.1.6, but getcwd error occurs.  Cannot install with Anaconda 4.3
    rooted_tree_pdf = directory + ".pdf"
    rooted_tree_svg = directory + ".svg"
    rooted_tree_path = "RAxML_bestTree.raxml.tre"
    if os.path.isfile("RAxML_bestTree.raxml.tre"):
        t = Tree(rooted_tree_path) #loads tree file
        ts = TreeStyle()
        for n in t.traverse():
            nstyle = NodeStyle()
            nstyle["size"] = 0 #removes dots from tree
            n.set_style(nstyle)
        def mylayout(node):
            if node.is_leaf():
                nameFace = AttrFace("name", fsize=9) #sets font size of leaves
                faces.add_face_to_node(nameFace, node, 0, position="branch-right")
        ts.layout_fn = mylayout #using custom layout above
        ts.show_leaf_name = False #using custom leaf size, so this is disabled
        ts.scale = 1000 #length of branches
        ts.branch_vertical_margin = 5 #spacing between branches
        ts.margin_left = 100
        ts.margin_right = 100
        ts.margin_top = 100
        ts.margin_bottom = 100
        t.render(rooted_tree_pdf, w=5000, tree_style=ts)
        t.render(directory + ".svg", w=500, tree_style=ts)
        os.rename(rooted_tree_path, directory + ".tre")

###################################################################
###################################################################
###################################################################

test_duplicate() #***FUNCTION CALL

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

elif options.species == "h37":
    
    qual_gatk_threshold = 150
    N_gatk_threshold = 200
    
    #Remove network path at and left of "Results"
    dependents_dir="/mycobacterium/tbc/h37/script_dependents/python/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/mycobacterium/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_000962.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/mycobacterium/tbc/h37/script2"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "suis1":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/suis1/script_dependents/python/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_017251-NC_017250.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/suis1/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "suis3":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/suis3/script_dependents/python/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NZ_CP007719-NZ_CP007718.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/suis3/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "suis4":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/suis4/script_dependents/python/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    #gbk_file = script_dependents + ""
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/suis4/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "ab1":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/abortus1/script_dependents/python/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_006932-NC_006933.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/abortus1/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "ab3":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/abortus3/script_dependents/python/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/CP007682-CP007683.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/abortus3/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "mel1":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/melitensis-bv1/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_003317-NC_003318.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/melitensis-bv1/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "mel2":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/melitensis-bv2/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_012441-NC_012442.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/melitensis-bv2/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "mel3":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/melitensis-bv3/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NZ_CP007760-NZ_CP007761.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/melitensis-bv3/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "canis":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/canis/script_dependents/python/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_010103-NC_010104.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/canis/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "ceti1":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/ceti1/script_dependents/python/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    #gbk_file = script_dependents + ""
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/ceti1/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif options.species == "ceti2":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/ceti2/script_dependents/python/script2"
    
    upload_to, remote, script_dependents = update_directory() #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_022905-NC_022906.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/ceti2/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if options.email == "s":
        email_list = "tod.p.stuber@usda.gov, jessica.a.hicks@usda.gov, christine.r.quance@usda.gov, suelee.robbe-austerman@aphis.usda.gov"

else:
    print ("EXIT AT SETTING OPTIONS")
    parser.print_help()
    sys.exit(0)

print ("\nSET VARIABLES")
print ("\tgenotypingcodes: %s " % genotypingcodes)
try:
    mygbk = True
    print ("\tgbk_file: %s " % gbk_file)
except NameError:
    mygbk = False
    print ("There is not a gbk file available")
print ("\tdefiningSNPs: %s " % definingSNPs)
print ("\texcelinfile: %s " % excelinfile)
print ("\tremove_from_analysis: %s " % remove_from_analysis)
print ("\tfilter_files: %s " % filter_files)
print ("\tbioinfoVCF: %s \n" % bioinfoVCF)
###

print ("\nChanging the VCF names")
names_not_changed = change_names() #***FUNCTION CALL

all_starting_files = glob.glob('starting_files/*vcf')
file_number = len(all_starting_files)

files = glob.glob('*vcf')
print ("REMOVING FROM ANALYSIS...")
wb = xlrd.open_workbook(remove_from_analysis)
ws = wb.sheet_by_index(0)
for each_sample in ws.col_values(0):
    each_sample = str(each_sample)
    each_sample = re.sub(r'(.*?)[._].*', r'\1', each_sample)
    myregex = re.compile(each_sample + '.*') # create regular expression to search for in VCF list
    for i in files:
        if myregex.search(i):
            print ("### --> %s removed from the analysis" % i)
            try:
                os.remove(i)
            except FileNotFoundError:
                print ("FileNotFoundError:")

#malformed_list = []
print ("CHECKING FOR EMPTY FILES...")
files = glob.glob('*vcf')
for i in files:
    if os.stat(i).st_size == 0:
        print ("### %s is an empty file and has been deleted" % i)
        error_list.append("File was empty %s" % i)
        os.remove(i)

print ("SORTING FILES...")
defining_snps = {}
inverted_position = {}
wb = xlrd.open_workbook(definingSNPs)
ws = wb.sheet_by_index(0)

for rownum in range(ws.nrows):
    position = ws.row_values(rownum)[1:][0]
    grouping = ws.row_values(rownum)[:1][0]
    # inverted positions will NOT be found in the passing positions
    # inverted positions are indicated in Defining SNPs by ending with "!"
    if position.endswith('!'):
        position = re.sub('!', '', position)
        inverted_position.update({position:grouping})
    else:
        defining_snps.update({position:grouping})
files = glob.glob('*vcf')

all_list_amb = {}
group_calls_list = []

print ("Grouping files...")
if options.debug_call:
    for i in files:
        dict_amb, group_calls, malformed = group_files(i)
        all_list_amb.update(dict_amb)
        group_calls_list.append(group_calls)
        error_list.append(malformed)
else:
    with futures.ProcessPoolExecutor() as pool:
        for dict_amb, group_calls, malformed in pool.map(group_files, files):
            all_list_amb.update(dict_amb)
            group_calls_list.append(group_calls) # make list of list
            error_list.append(malformed)
error_list = [x for x in error_list if x] # remove empty sets from list

print ("Getting directory list\n")
directory_list = next(os.walk('.'))[1] # get list of subdirectories
directory_list.remove('starting_files')

print ("Getting SNPs in each directory")
if options.debug_call:
    for i in directory_list:
        get_snps(i)
else:
    with futures.ProcessPoolExecutor(max_workers=cpu_count) as pool:
        pool.map(get_snps, directory_list)

runtime = (datetime.now() - startTime)
print ("\n\nruntime: %s:  \n" % runtime)

#############################################
#MAKE HTML FILE:
print ("<html>\n<head><style> table { font-family: arial, sans-serif; border-collapse: collapse; width: 40%; } td, th { border: 1px solid #dddddd; padding: 4px; text-align: left; font-size: 11px; } </style></head>\n<body style=\"font-size:12px;\">", file=htmlfile)
print ("<h2>Script ran using <u>%s</u> variables</h2>" % options.species.upper(), file=htmlfile)
print ("<h4>There are %s VCFs in this run</h4>" % file_number, file=htmlfile)

#OPTIONS
print ("Additional options ran: email: %s, inhouse: %s, options.filter: %s, all_vcf: %s, elite: %s, debug: %s, pilon: %s, uploaded: %s" % (options.email, options.inhouse, options.filter, options.all_vcf, options.elite, options.debug_call, options.pilon, options.upload), file=htmlfile)
if options.all_vcf:
    print ("\n<h4>All_VCFs is available</h4>", file=htmlfile)
elif options.elite:
    print ("\n<h4>Elite VCF comparison available</h4>", file=htmlfile)

#TIME
print ("\n<h4>Start time: %s <br>" % startTime, file=htmlfile)
print ("End time: %s <br>" % datetime.now(), file=htmlfile)
print ("Total run time: %s: </h4>" % runtime, file=htmlfile)

# ERROR LIST
if len(error_list) < 1:
    print ("<h2>No corrupt VCF removed</h2>", file=htmlfile)

else:
    print ("\n<h2>Corrupt VCF removed</h2>", file=htmlfile)
    for i in error_list:
        print ("%s <br>" % i, file=htmlfile)
    print ("<br>", file=htmlfile)

# AMBIGIOUS DEFINING SNPS
if len(all_list_amb) < 1:
    print ("\n<h2>No ambiguous defining SNPs</h2>", file=htmlfile)
else:
    print ("\n<h2>Defining SNPs are ambiguous.  They may be mixed isolates.</h2>", file=htmlfile)
    print ("<table>", file=htmlfile)
    print ("<tr align=\"left\"><th>Sample Name</th><th>Division</th><th>Absolute Position</th><tr>", file=htmlfile)
    ordered_all_list_amb = OrderedDict(sorted(all_list_amb.items()))
    for k, v in ordered_all_list_amb.items():
        k_split = k.split('\t')
        print ("<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (k_split[0], k_split[1], v), file=htmlfile)
    print ("</table>", file=htmlfile)
    print ("<br>", file=htmlfile)

#GROUPING TABLE
print ("<h2>Groupings</h2>", file=htmlfile)
print ("<table>", file=htmlfile)
print ("<tr align=\"left\"><th>Sample Name</th><tr>", file=htmlfile)
group_calls_list.sort(key=lambda x: x[0]) # sort list of list by first element
for i in group_calls_list:
    print ("<tr>", file=htmlfile)
    for x in i:
        print ("<td>%s</td>" % x, end='\t', file=htmlfile)
    print ("</tr>", file=htmlfile)
print ("</table>", file=htmlfile)

#FILES NOT RENAMED
if names_not_changed:
    print ("\n<h2>File names did not get changed:</h2>", file=htmlfile)
    for i in sorted(names_not_changed):
        print ("%s<br>" % i, file=htmlfile)

print ("</body>\n</html>", file=htmlfile)
#############################################
os.chdir(root)
zip("starting_files", "starting_files") # zip starting files directory
shutil.rmtree("starting_files")

htmlfile.close()

####send email:
def send_email():
    print ("Sending Email...")
    print ("Sending to:")

    msg = MIMEMultipart()
    msg['From'] = "tod.p.stuber@usda.gov"
    msg['To'] = email_list
    msg['Subject'] = "Script 2 " + options.species
    with open(htmlfile_name) as fp:
        msg.attach(MIMEText(fp.read(), 'html'))

    part = MIMEBase('application', "octet-stream")
    part.set_payload(open("summary_log.html", "r").read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', 'attachment; filename="summary_log.html"')
    msg.attach(part)

    smtp = smtplib.SMTP('10.10.8.12')
    smtp.send_message(msg)

    #smtp.send_message(msg)
    #smtp.send_message(msg.as_string())
    #smtp.sendmail(email_list, msg.as_string())
    #smtp.sendmail("tod.p.stuber@usda.gov", email_list, msg.as_string())
    smtp.quit()

if options.email == "none":
    print ("\n\temail not sent")
elif options.email:
    send_email()
    print ("\n\temail sent to: %s" % email_list)
else:
    print ("\n\temail not sent")

if options.upload:
    print ("Uploading Samples...")
    def copytree(src, dst, symlinks=False, ignore=None): #required to ignore permissions
        try:
            for item in os.listdir(src):
                s = os.path.join(src, item)
                d = os.path.join(dst, item)
                try:
                    if os.path.isdir(s):
                        shutil.copytree(s, d, symlinks, ignore)
                    else:
                        shutil.copy2(s, d)
                except shutil.Error:
                    pass
        except FileNotFoundError:
            print ("except FileNotFoundError: file not found")

    #upload to bioinfoVCF
    src = root
    dst = bioinfoVCF + "/" + os.path.basename(os.path.normpath(root))
    print ("\n\t%s is copying to %s" % (src, dst))
    os.makedirs(dst)
    copytree(src, dst)

print ("\n\tDONE\n")
#xvfb.kill

vdisplay.stop()

###############
# NOTES #######
###############
'''
Check python version by entering interpretter.
Install Anaconda if not already installed.
wget https://repo.continuum.io/archive/Anaconda3-4.0.0-Linux-x86_64.sh
bash Anaconda3-4.0.0-Linux-x86_64.sh

Script tested using Anaconda 4.0.0.  Other, newer, versions have shown not to work because some modules are incompatiable with Python >3.5.  If running a version other than 4.0.0 a new environment can be set without disrupting your current version (root) install.  Use Anaconda's default installation except when asked if to prepend to PATH choose yes.
If needing to create a new Anaconda environment with making changes to your current:
    $ conda create -n anaconda400 anaconda=4.0.0 anaconda
    To activate this environment, use:
    > source activate anaconda400
    To deactivate this environment, use:
    > source deactivate anaconda400
Log out and back in to your ssh session.
Go into the python interpretter by typing $ python
Your version will be Python 3.5.1 |Anaconda 4.0.0 (64-bit).  If it not, something did not run correctly.
Once Anaconda is installed
    $ easy_install ete3 
    $ easy_install pyvcf 
    $ easy_install biopython
    $ pip install xvfbwrapper
sudo yum install xorg-x11-server-Xvfb??? not sure this is needed

Helpful commands: conda env --help; conda env list
look at pandas version >>>pandas.__version__, cannot have pandas 0.19.2, pip install pandas==0.18.1, 4.0.0 should have 0.18.0

Run script and easy_install modules not available

X11:
    On Mac:
        Update with latest XQuartz??
        login with ssh -X
    On remote: 
        sudo yum install xorg-x11-server-Xvf python-xvfbwrapper xorg-x11-server-Xvfb.x86_6
        
MacOSX
download RAxML https://github.com/stamatak/standard-RAxML
Install in desired location
git clone https://github.com/stamatak/standard-RAxML.git
cd standard-RAxML/
make -f Makefile.SSE3.gcc
rm *.o
check PATH locations, $ echo $PATH, /usr/local/bin
ln -s /Users/Family/Documents/standard-RAxML/raxmlHPC-SSE3 ./raxml
check that it will be found in $PATH with $ which raxml

Linux
download RAxML https://github.com/stamatak/standard-RAxML
Install in desired location
git clone https://github.com/stamatak/standard-RAxML.git
cd standard-RAxML/
make -f Makefile.SSE3.gcc
rm *.o
check PATH locations, $ echo $PATH, /usr/local/bin
sudo ln -s /home/binfadmin/programs/standard-RAxML-master/raxmlHPC-SSE3 /usr/local/bin/raxml
check that it will be found in $PATH with $ which raxml



        
if bioinfo is not available cp script dependents to "/home/shared", or "/Users/Shared"
   for example, on a Mac, path to script dependents will look like: /Users/Shared/mycobacterium/tbc/tbbov/script_dependents/python/script2
 Use the tabs on the Excel Filter file to determine if chromosome name is correct
 Link from script1 into script2.  Chromosome name used in reference must be that used in script2.  See column one of a VCF.  This chromosome must be represent in the Defining SNPS, Filter file and gbk file.
 --> This shouldn't be needed if Anaconda 4.0 is being used: ete3 needs pyqt4, may need to first install pyqt4 and specify version using, $ conda install pyqt=4, ete3 will then install with easy_install
make sure xserver is running

# provide species when running script
Run: $ /path/to/script/script2 -s bovis

##Checklist
add a path setup for PC environment

'''
