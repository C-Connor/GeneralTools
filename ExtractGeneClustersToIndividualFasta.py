import os
import sys
from Bio import SeqIO
import glob
import argparse
import time

#parse command line arguments
cl_parser = argparse.ArgumentParser(description="Creates individual fasta files of gene sequences clustering file")
cl_parser.add_argument("-c","--clusters", help="Path to clustered_proteins file (produced by roary)",dest="clusters",type=str,required=True)
cl_parser.add_argument("-o","--output",help="Specify output path and directory name",type=str)
cl_parser.add_argument("-f","--fasta",help="Path to directories containing CDS fasta files. Please ensure all sequence files have same extension.",dest="fasta",type=str,required=True)
cl_parser.add_argument("-e","--extension",help="Specify fasta extenstion used for CDS files e.g. fasta, ffn, fa etc.",dest="extension",type=str,required=True)
cl_parser.add_argument("-l","--gene_list",help="File containing genes of interest, use gene names as they appear in the clustered_proteins files (same as roary pan genome matrix)",dest="genelist",type=str,required=True)

cl_args = cl_parser.parse_args()

#check input
if not os.path.isfile(cl_args.clusters):
    print('Could not find clustered_proteins file.')
    print('Exiting.')
    sys.exit()

if not os.path.isdir(cl_args.fasta):
    print('Could not find directory containing CDS fasta files.')
    print('Exiting.')
    sys.exit()

#if output is left empty make file name with time stamp
if cl_args.output == None:
    cl_args.output = os.path.join(os.getcwd(),'Clustered_Gene_Fasta_Files_'+time.strftime('%Y%m%d_%H%M')) + os.sep
elif not os.path.isdir(cl_args.output):
    os.makedirs(cl_args.output)

#get all ffn files from prokka directory
fastafiles = glob.glob(cl_args.fasta+'*.'+cl_args.extension)
if len(fastafiles) == 0:
    print('ERROR: Could not find any fasta files in specified directory!')
    print('Please check that the extension specified (.{0}) is correct.'.format(cl_args.extension))
    print('Exiting.')
    sys.exit()
    
print('Found {0} files matching extension: .{1}'.format(len(fastafiles),cl_args.extension))

#get group ids from gene list
groupidlist = dict()
with open(cl_args.genelist,'r') as genelist:
    for line in genelist:
        groupidlist[line.strip()]=''
        
#convert to set for faster searching
groupset = set(groupidlist.keys())
#use group ids to get individual gene ids from protein clustering file
with open(cl_args.clusters,'r') as clusterin:
    for line in clusterin:
        splitline = line.split(':')
        if splitline[0] in groupset:
            groupidlist[splitline[0]] = splitline[1].strip().split('\t')

print('{0} unique gene identifiers extracted.'.format(len(groupidlist.keys())))

#create index of fasta files
print('Indexing fasta files.')
fastaidx = SeqIO.index_db(':memory:', fastafiles, "fasta")

#write individual fasta files of nucleotide sequences
print('Writing gene sequences to fasta files.')
for gene in groupidlist.keys():
    SeqIO.write((fastaidx.get(sequence) for sequence in groupidlist.get(gene)),cl_args.output+str(gene).strip()+'.fasta','fasta')

print('Complete.')
