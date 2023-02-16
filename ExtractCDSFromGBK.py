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
    description="Extracts nucleotide or peptide sequences from CDS in GBK file. Can provide list of select genes to extract, provide EITHER gene name, locus tags or protein IDs.", 
    epilog="Example command:  ExtractCDSFromGBK.py nucleotide -i in.gbk -o out.fasta -l LocusTagList.txt"
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

group = cl_parser.add_mutually_exclusive_group()
group.add_argument("-g","--gene_list",
                       help="Filepath, list of CDS to extract, one per line, use GENE NAME e.g. 'thrL'. Default is None which extracts all CDS.",
                       dest='gene',type=str,default=None,metavar="FILEPATH")
group.add_argument("-l","--locus_list",
                       help="Filepath, list of CDS to extract, one per line, use LOCUS TAG e.g. 'tag_001'. Default is None which extracts all CDS.",
                       dest='locus_tag',type=str,default=None,metavar="FILEPATH")
group.add_argument("-p","--protid_list",
                       help="Filepath, list of CDS to extract, one per line, use PROTEIN ID e.g. 'AB1234.1'. Default is None which extracts all CDS.",
                       dest='protein_id',type=str,default=None,metavar="FILEPATH")


cl_args = cl_parser.parse_args()

#check user output file exists or make directories if necessary
if cl_args.output.exists():
    raise FileExistsError(f"Output file already exists: {cl_args.output}")
if not cl_args.output.parent.exists():
    cl_args.output.parent.mkdir(parents=True,exist_ok=False)


#check if CDS list has been provided and which format
if any([vars(cl_args)[v] != None for v in ['gene','locus_tag','protein_id']]):
    #cds set
    for arg,val in vars(cl_args).items():
        if val != None:
            cds = (arg,val) #relies on mutually exclusive options
else:
    cds= None

#read in selected CDS if set
if cds != None:
    with open(cds[1]) as infile:
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
    if cds != None:
        name = f.qualifiers[cds[0]][0]
    else:
        name = f.qualifiers['protein_id'][0]
    return SeqRecord(
        seq,
        id="|".join([y for x in [f.qualifiers[key] for key in f.qualifiers.keys() if key in set(['protein_id','gene','db_xref'])] for y in x]),
        name=name, 
        description = f.qualifiers['product'][0]
    )

#iterate over CDS features in GBK file calling Pull function 
seqs = []
for f in genome.features:
    if (f.type == 'CDS') and ('protein_id' in f.qualifiers.keys()):
        if cds == None:
            seqs.append(Pull(f))
        else:
            if (cds[0] in f.qualifiers.keys()) and (f.qualifiers[cds[0]][0] in select_cds):
                seqs.append(Pull(f))

#check if all CDS have been found, list those that haven't
if cds != None:
    if len(seqs) == 0:
        print('ERROR: None of the specified CDS not found in GBK file.\nExiting.')
        sys.exit(1)
    elif len(seqs) < len(select_cds):
        print(f"Not all Selected CDS specified in '{cds[1]}' were found in '{cl_args.gbkin}'")
        print(f"Missing CDS: {set(select_cds) - set(x.name for x in seqs)}")

#write found sequences to file in fasta format
with open(cl_args.output,'w') as outfile:
    SeqIO.write(seqs,outfile,'fasta')
