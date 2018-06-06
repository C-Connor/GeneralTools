#!/usr/bin/env python2
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import sys
import tempfile
import argparse
import time

#parse command line arguments
cl_parser = argparse.ArgumentParser(description="Translates nucelotide fasta sequences into amino acid multifasta using bacterial genetic table.")
cl_parser.add_argument("-i","--fasta_in", help="Fasta file containing nucleotide sequences to be translated",dest="fasta_in",type=str,required=True)
cl_parser.add_argument("-o","--output",help="Specify output path and file name",dest="output",type=str)
cl_args = cl_parser.parse_args()

#function for translating nucleotide to amino acid
def make_protein_record(nuc_record):
	try:
		return SeqRecord(seq = nuc_record.seq.translate(cds=True,table='Bacterial'), id = nuc_record.id, description = '')
	except Exception as ex:
		if type(ex).__name__ == 'TranslationError':
			return ex
		else:
			raise

#create function to give to SeqIO for creating dictionary keys from pangenome fasta file
def inputfastakey(identifier):
    part = identifier.split('|')[1]
    return(part)

if not os.path.isfile(cl_args.fasta_in):
    print('Could not find input fasta file')
    print('Exiting')
    sys.exit()

if cl_args.output != None:
    try:
        os.makedirs(os.path.abspath(os.path.dirname(cl_args.output)))
    except OSError as er:
        if er.errno == os.errno.EEXIST:
            print('Output directory already exists.')
        else:
            raise
if cl_args.output == None:
    cl_args.output = os.getcwd()+os.sep+'CDStoAA_'+time.strftime('%Y%m%d_%H%M')+'.fasta'

#if using roary pan genome fasta, modify fasta file to use | symbols in comment lines
fastatemp = tempfile.NamedTemporaryFile(suffix = '.fasta', prefix = 'cds_temp', dir = os.path.dirname(cl_args.output), bufsize=0)
with open(cl_args.fasta_in,'r') as rawfasta:
    rawfasta.seek(0)
    for line in rawfasta:
        newline = line.replace(' ','|')
        fastatemp.write(newline)

print('Translating CDS using bacterial table.')
proteins = []
errorlog = []
errorcount = 1
for nuc_rec in SeqIO.parse(fastatemp.name,"fasta"):
	peptide = make_protein_record(nuc_rec)
	if type(peptide) is not SeqRecord:
		errorlog.append('Error {0}:	Nucleotide ID: {1}	Error Message: {2}\n'.format(errorcount,nuc_rec.id,peptide))
		errorcount += 1
	else:
		proteins.append(peptide)


print('Writing sequences to file.')
SeqIO.write(proteins, cl_args.output, "fasta")

fastatemp.close()

print('Sequences written to {0}'.format(cl_args.output))

if len(errorlog) > 0:
	print('WARNING: Translation errors were encountered. See log file for details.')
	print('Logfile location: '+os.path.dirname(cl_args.output)+'CDSToAminoAcid_logfile_'+time.strftime('%Y%m%d_%H%M')+'.txt')
	with open(os.path.join(os.path.dirname(cl_args.output),'CDSToAminoAcid_logfile_'+time.strftime('%Y%m%d_%H%M')+'.txt'),'w') as logfile:
		for entry in errorlog:
			logfile.write(str(entry))
