#!/usr/bin/env python2
from Bio import SeqIO
import argparse
import os
import sys
import time

#command line arguments
cl_parser = argparse.ArgumentParser(description="Extracts locus tag, gene name, gene product, protein accession numbers, start and end coordinates of features from genbank files.")
cl_parser.add_argument("-i","--gbk_in", help="Genbank file to extract CDS features from",dest="gbkin",type=str,required=True)
cl_parser.add_argument("-o","--output",help="Name and path to put output file",dest="output",type=str)
cl_parser.add_argument("-l","--list_features",help="Print list of features which will be extracted from genbank file and exit",dest='flist',action='store_true')
cl_args = cl_parser.parse_args()

feature_list = ['CDS','tRNA','rRNA','misc_RNA','mobile_element','ncRNA','repeat_region']


#check user input
if cl_args.flist:
	print('Features to be extracted from genbank file:')
	print('')
	for i in feature_list:
		print(i)
	print('')
	print('Exiting')
	sys.exit()

if not os.path.isfile(cl_args.gbkin):
    print('Could not find input genbank file')
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
    cl_args.output = os.getcwd()+os.sep+'FeaturesFromGBK_'+time.strftime('%Y%m%d_%H%M')+'.tab'

#extract CDS features and info, and write to file
print('Extracting  features...')
with open(cl_args.output,'w') as outhandle:
    outhandle.write('Feature_type\tLocus_tag\tgene_name\tproduct\tdb_xref\tprotein_id\tstart\tend\n')
    for record in SeqIO.parse(cl_args.gbkin,'genbank'):
        for seqrecord in record.features:
			if seqrecord.type in set(feature_list):
				try:
					locus_tag = seqrecord.qualifiers.get('locus_tag')[0]
				except TypeError:
					locus_tag = seqrecord.qualifiers.get('locus_tag')
				try:
                			gene = seqrecord.qualifiers.get('gene')[0]
				except TypeError:
					gene = seqrecord.qualifiers.get('gene')
				try:
					product = seqrecord.qualifiers.get('product')[0]
				except TypeError:
					product = seqrecord.qualifiers.get('product')
				try:
					dbxref = seqrecord.qualifiers.get('db_xref')[0]
				except TypeError:
					dbxref = seqrecord.qualifiers.get('db_xref')
				try:
					proteinid = seqrecord.qualifiers.get('protein_id')[0]
				except TypeError:
					proteinid = seqrecord.qualifiers.get('protein_id')
				start = str(seqrecord.location.start)
				end = str(seqrecord.location.end)
				feature_type = seqrecord.type
				outhandle.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(feature_type,locus_tag,gene,product,dbxref,proteinid,start,end))

print('Complete')
print('Output file in {0}'.format(cl_args.output))
