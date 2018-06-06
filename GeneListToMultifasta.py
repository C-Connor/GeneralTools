#!/usr/bin/env python2
import os
import sys
from Bio import SeqIO
import tempfile
import time
import argparse

#parse command line arguments
cl_parser = argparse.ArgumentParser(description="Extract sequences of interest from multifasta file and write to new multifasta file.")
cl_parser.add_argument("-f","--fasta_in", help="Input multifasta file",dest="fastain",type=str,required=True)
cl_parser.add_argument("-o","--output",help="Specify output path",dest="output",type=str)
cl_parser.add_argument("-g","--gene_list", help="List of genes of interest, one per line.",dest="genelist",type=str,required=True)

cl_args = cl_parser.parse_args()

#check user input
if not os.path.isfile(cl_args.fastain):
    print('Cannot find input fasta file.')
    print('Exiting')

if not os.path.isfile(cl_args.genelist):
    print('Cannot find gene list')
    print('Exiting')

if cl_args.output != None:
    try:
        os.makedirs(os.path.abspath(os.path.dirname(cl_args.output)))
    except OSError as er:
        if er.errno == os.errno.EEXIST:
            print('Output directory already exists.')
        else:
            raise
if cl_args.output == None:
    cl_args.output = os.getcwd()+os.sep+'Multifasta_'+time.strftime('%Y%m%d_%H%M')+'.fasta'


#create function to give to SeqIO for creating dictionary keys from pangenome fasta file
def inputfastakey(identifier):
    part = identifier.split('|')[1]
    return(part)

genelist = list()
with open(cl_args.genelist,'r') as ingenelist:
    for line in ingenelist:
        genelist.append(line.strip())

#check if fasta file is formatted correctly?

#create temporary file to store modified reference fasta
#modify reference fasta to identifier line contains | between genome name and gene name
fastatemp = tempfile.NamedTemporaryFile(suffix = '.fasta', prefix = 'pan_gen_temp', dir = os.path.dirname(cl_args.output), bufsize=0)
with open(cl_args.fastain,'r') as rawfasta:
    rawfasta.seek(0)
    for line in rawfasta:
        newline = line.replace(' ','|')
        fastatemp.write(newline)

#create index of temporary fasta file using gene name as dictionary key
reference_fasta_dict = SeqIO.index(fastatemp.name,'fasta', key_function=inputfastakey)

#write fasta sequences to output files
#create generator / iterator for indexed fasta file using gene name lists, then pass to SeqIO output
try:
    FastaRecords = (reference_fasta_dict[gene] for gene in genelist)
    SeqIO.write(FastaRecords, cl_args.output, 'fasta')
except KeyError as err:
    #print(err)
    print('-------------')
    print('KeyError: {0}'.format(gene))
    print('Input gene list contains an entry not present in the provided fasta file')
    print('List of genes in file must be separated by newline "\\n" characters only with no extra blank lines.')
    print('-------------')
    print('Output file is incomplete.')
    print('Exiting')
    fastatemp.close()
    reference_fasta_dict.close()
    sys.exit()

#close file connections
#closing temporary file will delete file
reference_fasta_dict.close()
fastatemp.close()
print('Output in: {0}'.format(cl_args.output))
print('Complete')
