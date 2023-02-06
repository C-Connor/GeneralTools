#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import os
import errno
import sys
import time
import pathlib


#command line arguments
cl_parser = argparse.ArgumentParser(
    description="Extracts nucleotide or peptide sequences from CDS in GBK file.", 
    epilog="Example command:  ExtractCDSFromGBK.py nucleotide -i in.gbk -o out.fasta -c CDSlist.txt"
    )
cl_parser.add_argument("-i","--gbk_in", 
                       help="Genbank file to extract CDS sequences from",
                       dest="gbkin",type=pathlib.Path,required=True)
cl_parser.add_argument("-o","--output",
                       help="Name and path to put output file",
                       dest="output",type=pathlib.Path)
cl_parser.add_argument("SeqType",
                       help="Extract nucleotide or peptide sequences",
                       choices=['nucleotide','peptide'])
cl_parser.add_argument("-c","--cds_list",
                       help="Filepath, list of CDS to extract, one per line, use gene name e.g. thrL. Default is None which extracts all CDS.",
                       dest='cds',type=str,default=None)
cl_args = cl_parser.parse_args()

#check user output file exists or make directories if necessary
if cl_args.output.exists():
    raise FileExistsError(f"Output file already exists: {cl_args.output}")
if not cl_args.output.parent.exists():
    cl_args.output.parent.mkdir(parents=True,exist_ok=False)

#read in selected CDS if set
if cl_args.cds != None:
    with open(cl_args.cds) as infile:
        select_cds = [line.strip() for line in infile]

#read GBK file
with open(cl_args.gbkin,'r') as infile:
    genome = next(SeqIO.parse(infile,'genbank')) #not sure about this, what about multiple records/genomes in genbank?

#function to pull sequence info from GBK file
def Pull(f):
    if cl_args.SeqType == 'peptide':
        seq = Seq(f.qualifiers['translation'][0])
    else:
        seq = f.extract(genome.seq)
    return SeqRecord(
        seq,
        id="|".join([f.qualifiers['protein_id'][0],f.qualifiers['gene'][0]]+f.qualifiers['db_xref']),
        name=f.qualifiers['gene'][0],
        description = f.qualifiers['product'][0]
    )

#iterate over CDS features in GBK file calling Pull function 
seqs = []
for f in genome.features:
    if (f.type == 'CDS') and ('protein_id' in f.qualifiers.keys()):
        if cl_args.cds == None:
            seqs.append(Pull(f))
        else:
            if f.qualifiers['gene'][0] in select_cds:
                seqs.append(Pull(f))

#check if all CDS have been found, list those that haven't
if len(seqs) == 0:
    print('ERROR: Specified CDS not found in GBK file')
    sys.exit(1)
elif len(seqs) < len(select_cds):
    print(f"Not all Selected CDS specified in '{cl_args.cds}' were found in '{cl_args.gbkin}'")
    print(f"CDS not found: {set(select_cds) - set(x.name for x in seqs)}")

#write found sequences to file in fasta format
with open(cl_args.output,'w') as outfile:
    SeqIO.write(seqs,outfile,'fasta')
