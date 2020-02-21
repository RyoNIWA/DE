#! /usr/bin/env python3

import subprocess
import glob
from optparse import OptionParser
import os
import shutil
import sys

class pycolor:
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    END = '\033[0m'

parser=OptionParser(description='This script was made by Ryo NIWA from Gifu University. Thank you for using DifferenceExtaction. After finishing analysis, please put the unique genes to BLASTkoala to utilize the data!')
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~                  ~~~")
X = pycolor.BLUE + " Created by Ryo " + pycolor.END
Y = "~~~ {} ~~~"
Z = Y.format(X)
print(Z)
print("~~~  Gifu University ~~~")
print("~~~                  ~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print()
parser.add_option('-1', dest="input1", help='DifferenceExtraction accepts 2 files of faa. Only Amino acid sequences are applicable currently. The unique gene sets of input1 will be stored in INPUT1 dirrectory.')
parser.add_option('-2', dest="input2", help='DifferenceExtraction accepts 2 files of faa. Only Amino acid sequences are applicable currently. The unique gene sets of input2 will be stored in INPUT2 dirrectory.')
parser.add_option('-e', dest="evalue", help='Evalue may change the results of BLAST. You can define the threshold. The number of each core genes and unique genes will be differred depending on e-value of BLAST. default=0.0001', default=0.0001)
(options, args) = parser.parse_args()
print("This script was made by Ryo NIWA from Gifu University. Thank you for using DifferenceExtaction. After finishing analysis, please put the unique genes to BLASTkoala to utilize the data!")
print()
print()
#PATH check
try:
    print("Check the path to input files:")
    input1_faa = os.path.abspath(options.input1)
    input2_faa = os.path.abspath(options.input2)
    print(input1_faa)
    print(input2_faa)
    print()
except:
    Error1 = pycolor.RED + "Please define input1 and input2" + pycolor.END
    print(Error1)
    sys.exit()

try:
    print("Check the path to BLAST:")
    print(shutil.which('makeblastdb'))
    print()
except:
    Error1 = pycolor.RED + "Please install BLAST using brew or conda." + pycolor.END
    print(Error1)
    sys.exit()

OK = pycolor.BLUE + "All DEPENDENCIES are OK!" + pycolor.END
print(OK)
print()
print()

# This process makes DATABASE directory
new_dir_path_recursive1 = './BLAST/DATABASE'
os.makedirs(new_dir_path_recursive1, exist_ok=True)

# This process make Database file for BLAST
First = pycolor.GREEN + "Searching for unique genes in input1 will start soon!" + pycolor.END
print(First)
print("Let us make database for BLASTp")
try:
    makeblastdb = shutil.which('makeblastdb')
    blastdb_cmd1_A = 'makeblastdb -in {} -dbtype prot -out ./BLAST/DATABASE/input1 -parse_seqids -logfile ./BLAST/DATABASE/makeblastdb2.log'
    blastdb_cmd1 = blastdb_cmd1_A.format(input1_faa)
    DB_process = subprocess.Popen(blastdb_cmd1,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    print("Operation starts now!")
    print()
    print()
    blastdb_cmd2_A = 'makeblastdb -in {} -dbtype prot -out ./BLAST/DATABASE/input2 -parse_seqids -logfile ./BLAST/DATABASE/makeblastdb1.log'
    blastdb_cmd2 = blastdb_cmd2_A.format(input2_faa)
    DB_process = subprocess.Popen(blastdb_cmd2,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    DB_process.wait()
except:
    Error2 = pycolor.RED + "makeblastdb was failed to complete. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error2)
    sys.exit()

# # This process makes result_input1 directory
new_dir_path_recursive2 = './BLAST/result_input1'
os.makedirs(new_dir_path_recursive2, exist_ok=True)

#blastp
print("BLASTp may take long time to complete.")
Query = input1_faa
Database = input2_faa
info = "1st BLASTp [Query: {}, Database: {}]"
infom = info.format(Query, Database)
print(infom)
try:
    blastp = shutil.which('blastp')
    blastp_cmd_A = 'blastp -query {} -db ./BLAST/DATABASE/input2 -evalue {} -outfmt {} > ./BLAST/result_input1/blastp.outfmt7'
    evalue = options.evalue
    blast_set = '"7 qseqid"'
    blastp_cmd = blastp_cmd_A.format(input1_faa, evalue, blast_set)
    BLASTP = subprocess.Popen(blastp_cmd,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    BLASTP.wait()
except:
    Error3 = pycolor.RED + "BLASTP was failed to complete. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error3)
    sys.exit()

#Extraction of ID of each alignment
print("Complete 1st BLASTp, and now move on to finding differences between query and subject")
print("Extracting the best hits from the result of blastp")
print()
try:
    extract_A = 'cat ./BLAST/result_input1/blastp.outfmt7 | awk {} | grep -v {} > ./BLAST/result_input1/best_hits.txt'
    awk = '"/hits found/{getline;print}¥"'
    awk_set ='"#"'
    extract = extract_A.format(awk, awk_set)
    res = subprocess.Popen(extract,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    res.wait()
except:
    Error4 = pycolor.RED + "Error was occured while you extract IDs of each alignment. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error4)
    sys.exit()

# Divide core genes and unique genes
# core
print()
print("Divide core genes and unique genes")
try:
    devision_A = 'grep _ {} > ./BLAST/result_input1/list.txt'
    devision = devision_A.format(input1_faa)
    devision_cmd = subprocess.Popen(devision,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    devision_cmd.wait()
except:
    Error5 = pycolor.RED + "Error was occured while you divide core genes and unique genes of core genes. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error5)
    sys.exit()

try:
    devision_A = 'cat ./BLAST/result_input1/list.txt | grep -f ./BLAST/result_input1/best_hits.txt > ./BLAST/result_input1/core_genes_trash.txt'
    devision = devision_A.format(input1_faa)
    devision_cmd = subprocess.Popen(devision,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    devision_cmd.wait()
except:
    Error6 = pycolor.RED + "Error was occured while you perform command 'cat and grep'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error6)
    sys.exit()

try:
    sed = 'sed -e {} ./BLAST/result_input1/core_genes_trash.txt > ./BLAST/result_input1/core_genes.txt'
    sed_set = "'s/>//g'"
    sed_sed = sed.format(sed_set)
    sed_cmd = subprocess.Popen(sed_sed,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    sed_cmd.wait()
except:
    Error7 = pycolor.RED + "Error was occured while you perform command 'sed'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error7)
    sys.exit()

print("Core genes were extracted from input1.faa successfully.")
print()
# unique
try:
    devision_A = 'grep _ {} > ./BLAST/result_input1/list2.txt'
    devision = devision_A.format(input2_faa)
    devision_cmd = subprocess.Popen(devision,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    devision_cmd.wait()
except:
    Error8 = pycolor.RED + "Error was occured while you divide core genes and unique genes of unique genes. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error8)
    sys.exit()

try:
    devision_A = 'cat ./BLAST/result_input1/list.txt | grep -v -f ./BLAST/result_input1/best_hits.txt > ./BLAST/result_input1/unique_genes_trash.txt'
    devision = devision_A.format(input1_faa)
    devision_cmd = subprocess.Popen(devision,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    devision_cmd.wait()
except:
    Error9 = pycolor.RED + "Error was occured while you perform command 'cat and grep'. Please check whether you satisfy all dependencies" + pycolor.END
    print(Error9)
    sys.exit()

try:
    sed = 'sed -e {} ./BLAST/result_input1/unique_genes_trash.txt > ./BLAST/result_input1/unique_genes.txt'
    sed_set = "'s/>//g'"
    sed_sed = sed.format(sed_set)
    sed_cmd = subprocess.Popen(sed_sed,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    sed_cmd.wait()
except:
    Error10 = pycolor.RED + "Error was occured while you perform command 'sed'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error10)
    sys.exit()

print("Unique genes were extracted from input1.faa successfully.")
print("Final Step is installing sequences of amino acid sequences as faa file in the folder of AminoAcidSequence_input1")

#Remove non-necessary files
os.remove('./BLAST/result_input1/unique_genes_trash.txt')
os.remove('./BLAST/result_input1/core_genes_trash.txt')

print("Preparation of installation of sequences of amino acid is done. Operation will be completed soon.")
aminoacid1 = './AminoAcidSequence_input1'
os.makedirs(aminoacid1, exist_ok=True)

try:
    cat = 'cat ./BLAST/result_input1/core_genes.txt | cut -d " " -f 1 > ./AminoAcidSequence_input1/core_genesID.txt'
    cat_cmd = subprocess.Popen(cat,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    cat_cmd.wait()
except:
    Error11 = pycolor.RED + "Error was occured while you perform command 'cat'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error11)
    sys.exit()

try:
    cat = 'cat ./BLAST/result_input1/unique_genes.txt | cut -d " " -f 1 > ./AminoAcidSequence_input1/unique_genesID.txt'
    cat_set = '" "'
    cat_cat = cat.format(cat_set)
    cat_cmd = subprocess.Popen(cat_cat,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    cat_cmd.wait()
except:
    Error12 = pycolor.RED + "Error was occured while you perform command 'cat'. Please check whether you satisfy all dependencies. 2" + pycolor.END
    print(Error12)
    sys.exit()

try:
    blastdbcmd = shutil.which('blastdbcmd')
    blastdb = 'blastdbcmd -entry_batch ./AminoAcidSequence_input1/core_genesID.txt -db ./BLAST/DATABASE/input1 -dbtype prot -out ./AminoAcidSequence_input1/core_genes.faa'
    blastdbcmd_cmd = subprocess.Popen(blastdb,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    blastdbcmd_cmd.wait()
except:
    Error13 = pycolor.RED + "Error was occured while you perform command 'blastdbcmd'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error13)
    sys.exit()

try:
    blastdbcmd = shutil.which('blastdbcmd')
    blastdbcmd2 = 'blastdbcmd -entry_batch ./AminoAcidSequence_input1/unique_genesID.txt -db ./BLAST/DATABASE/input1 -dbtype prot -out ./AminoAcidSequence_input1/unique_genes.faa'
    blastdbcmd2_cmd = subprocess.Popen(blastdbcmd2,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    blastdbcmd2_cmd.wait()
except:
    Error14 = pycolor.RED + "Error was occured while you perform command 'blastdbcmd'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error14)
    sys.exit()

# 2nd step
print()
print()
Half = pycolor.GREEN + "A half of the whole analysis was done! Wait a little bit more till the end!" + pycolor.END
print(Half)
Half2 = pycolor.GREEN + "Searching for unique genes in input2 will start soon!" + pycolor.END
print(Half2)
new_dir_path_recursive3 = './BLAST/result_input2'
os.makedirs(new_dir_path_recursive3, exist_ok=True)
print("BLASTp will be performed again.")
input2_faa = os.path.abspath(options.input1)
input1_faa = os.path.abspath(options.input2)
Query = input1_faa
Database = input2_faa
info = "2nd BLASTp [Query: {}, Database: {}]"
infom = info.format(Query, Database)
print(infom)

try:
    makeblastdb = shutil.which('makeblastdb')
    blastdb_cmd1_A = 'makeblastdb -in {} -dbtype prot -out ./BLAST/DATABASE/input1 -parse_seqids -logfile ./BLAST/DATABASE/makeblastdb2.log'
    blastdb_cmd1 = blastdb_cmd1_A.format(input1_faa)
    DB_process = subprocess.Popen(blastdb_cmd1,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    print("Operation starts now!")
    print()
    print()
    blastdb_cmd2_A = 'makeblastdb -in {} -dbtype prot -out ./BLAST/DATABASE/input2 -parse_seqids -logfile ./BLAST/DATABASE/makeblastdb1.log'
    blastdb_cmd2 = blastdb_cmd2_A.format(input2_faa)
    DB_process = subprocess.Popen(blastdb_cmd2,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    DB_process.wait()
except:
    Error2 = pycolor.RED + "makeblastdb was failed to complete. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error2)
    sys.exit()

try:
    blastp = shutil.which('blastp')
    blastp_cmd_A = 'blastp -query {} -db ./BLAST/DATABASE/input2 -evalue {} -outfmt {} > ./BLAST/result_input2/blastp.outfmt7'
    evalue = options.evalue
    blast_set = '"7 qseqid"'
    blastp_cmd = blastp_cmd_A.format(input1_faa, evalue, blast_set)
    BLASTP = subprocess.Popen(blastp_cmd,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    BLASTP.wait()
except:
    Error3 = pycolor.RED + "BLASTP was failed to complete. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error3)
    sys.exit()

#Extraction of ID of each alignment
print("Complete 1st BLASTp, and now move on to finding differences between query and subject")
print("Extracting the best hits from the result of blastp")
print()
try:
    extract_A = 'cat ./BLAST/result_input2/blastp.outfmt7 | awk {} | grep -v {} > ./BLAST/result_input2/best_hits.txt'
    awk = '"/hits found/{getline;print}¥"'
    awk_set ='"#"'
    extract = extract_A.format(awk, awk_set)
    res = subprocess.Popen(extract,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    res.wait()
except:
    Error4 = pycolor.RED + "Error was occured while you extract IDs of each alignment. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error4)
    sys.exit()

# Divide core genes and unique genes
# core
print()
print("Divide core genes and unique genes")
try:
    devision_A = 'grep _ {} > ./BLAST/result_input2/list.txt'
    devision = devision_A.format(input1_faa)
    devision_cmd = subprocess.Popen(devision,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    devision_cmd.wait()
except:
    Error5 = pycolor.RED + "Error was occured while you divide core genes and unique genes of core genes. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error5)
    sys.exit()

try:
    devision_A = 'cat ./BLAST/result_input2/list.txt | grep -f ./BLAST/result_input2/best_hits.txt > ./BLAST/result_input2/core_genes_trash.txt'
    devision = devision_A.format(input1_faa)
    devision_cmd = subprocess.Popen(devision,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    devision_cmd.wait()
except:
    Error6 = pycolor.RED + "Error was occured while you perform command 'cat and grep'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error6)
    sys.exit()

try:
    sed = 'sed -e {} ./BLAST/result_input2/core_genes_trash.txt > ./BLAST/result_input2/core_genes.txt'
    sed_set = "'s/>//g'"
    sed_sed = sed.format(sed_set)
    sed_cmd = subprocess.Popen(sed_sed,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    sed_cmd.wait()
except:
    Error7 = pycolor.RED + "Error was occured while you perform command 'sed'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error7)
    sys.exit()

print("Core genes were extracted from input2.faa successfully.")
print()
# unique
try:
    devision_A = 'grep _ {} > ./BLAST/result_input2/list2.txt'
    devision = devision_A.format(input2_faa)
    devision_cmd = subprocess.Popen(devision,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    devision_cmd.wait()
except:
    Error8 = pycolor.RED + "Error was occured while you divide core genes and unique genes of unique genes. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error8)
    sys.exit()

try:
    devision_A = 'cat ./BLAST/result_input2/list.txt | grep -v -f ./BLAST/result_input2/best_hits.txt > ./BLAST/result_input2/unique_genes_trash.txt'
    devision = devision_A.format(input1_faa)
    devision_cmd = subprocess.Popen(devision,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    devision_cmd.wait()
except:
    Error9 = pycolor.RED + "Error was occured while you perform command 'cat and grep'. Please check whether you satisfy all dependencies" + pycolor.END
    print(Error9)
    sys.exit()

try:
    sed = 'sed -e {} ./BLAST/result_input2/unique_genes_trash.txt > ./BLAST/result_input2/unique_genes.txt'
    sed_set = "'s/>//g'"
    sed_sed = sed.format(sed_set)
    sed_cmd = subprocess.Popen(sed_sed,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    sed_cmd.wait()
except:
    Error10 = pycolor.RED + "Error was occured while you perform command 'sed'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error10)
    sys.exit()

print("Unique genes were extracted from input2.faa successfully.")
print("Final Step is installing sequences of amino acid sequences as faa file in the folder of AminoAcidSequence_input2")

#Remove non-necessary files
os.remove('./BLAST/result_input2/unique_genes_trash.txt')
os.remove('./BLAST/result_input2/core_genes_trash.txt')

print("Preparation of installation of sequences of amino acid is done. Operation will be completed soon.")
aminoacid2 = './AminoAcidSequence_input2'
os.makedirs(aminoacid2, exist_ok=True)

try:
    cat = 'cat ./BLAST/result_input2/core_genes.txt | cut -d " " -f 1 > ./AminoAcidSequence_input2/core_genesID.txt'
    cat_cmd = subprocess.Popen(cat,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    cat_cmd.wait()
except:
    Error11 = pycolor.RED + "Error was occured while you perform command 'cat'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error11)
    sys.exit()

try:
    cat = 'cat ./BLAST/result_input2/unique_genes.txt | cut -d " " -f 1 > ./AminoAcidSequence_input2/unique_genesID.txt'
    cat_set = '" "'
    cat_cat = cat.format(cat_set)
    cat_cmd = subprocess.Popen(cat_cat,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    cat_cmd.wait()
except:
    Error12 = pycolor.RED + "Error was occured while you perform command 'cat'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error12)
    sys.exit()

try:
    blastdbcmd = shutil.which('blastdbcmd')
    blastdb = 'blastdbcmd -entry_batch ./AminoAcidSequence_input2/core_genesID.txt -db ./BLAST/DATABASE/input1 -dbtype prot -out ./AminoAcidSequence_input2/core_genes.faa'
    blastdbcmd_cmd = subprocess.Popen(blastdb,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    blastdbcmd_cmd.wait()
except:
    Error13 = pycolor.RED + "Error was occured while you perform command 'blastdbcmd'. Please check whether you satisfy all dependencies." + pycolor.END
    print(Error13)
    sys.exit()

try:
    blastdbcmd = shutil.which('blastdbcmd')
    blastdbcmd2 = 'blastdbcmd -entry_batch ./AminoAcidSequence_input2/unique_genesID.txt -db ./BLAST/DATABASE/input1 -dbtype prot -out ./AminoAcidSequence_input2/unique_genes.faa'
    blastdbcmd2_cmd = subprocess.Popen(blastdbcmd2,
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
    blastdbcmd2_cmd.wait()
except:
    Error14 = pycolor.RED + "Error was occured while you perform command 'blastdbcmd'. Please check whether you satisfy all dependencies. 3" + pycolor.END
    print(Error14)
    sys.exit()

print()
print()
Final = pycolor.GREEN + "Thank you for using DifferenceExtraction. Your analysis was successfully done!" + pycolor.END
print(Final)
print()
print()

import sys
line = 0
# ファイルを開く
with open('./AminoAcidSequence_input2/core_genesID.txt', "r") as f:
    # 一行ずつ読み込む
    for data in f:
        # 行数を加算
        line += 1
line1 = 0
# ファイルを開く
with open('./AminoAcidSequence_input2/unique_genesID.txt', "r") as f:
    # 一行ずつ読み込む
    for data in f:
        # 行数を加算
        line1 += 1
line2 = 0
# ファイルを開く
with open('./AminoAcidSequence_input1/unique_genesID.txt', "r") as f:
    # 一行ずつ読み込む
    for data in f:
        # 行数を加算
        line2 += 1
        
Final = pycolor.YELLOW + "Score Reports" + pycolor.END
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(Final)
print()
print("The number of genes in input1: " + str(line+line2))
print("The number of genes in input2: " + str(line+line1))
print("Core Genes: " + str(line))
print("Unique genes in input1:" + str(line2) +  " Percentage: " + str(line2/(line+line2)*100) + "%")
print("Unique genes in input2:" + str(line1) +  " Percentage: " +  str(line1/(line+line1)*100) + "%")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print()
print()
print("The definition of Percentage is [Number of unique genes of input1 or input2]/[Number of total genes of input1 or input2]")
print()
print()
print("Operation  completed!")
print()
