
"""
The goal for this program is to provide an iterable function that can be run
across a reference genome to isolate and return sequences that lie between two
primer sequences provided, given certain variables.
The function takes the following files:

    The reference genome path string 
        Either as a fasta or fastq file ending in that file name

    The primers file that contains both the forward and reverse primer sequences

    Parameters for determining whether the region selected is ITS or TEF
        Implemented as a maximum and minimum length of the extracted sequence
            This allows for extraction of other regions

"""

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import subprocess
import os
import argparse
import fileinput

parser = argparse.ArgumentParser(description="""
The goal for this program is to provide an iterable function that can be run
across a reference genome to isolate and return sequences that lie between two
primer sequences provided, given certain variables.
The function takes the following files:

    The reference genome path string 
        Either as a fasta or fastq file ending in that file name

    The primers file that contains both the forward and reverse primer sequences

    Parameters for determining whether the region selected is ITS or TEF
        Implemented as a maximum and minimum length of the extracted sequence
            This allows for extraction of other regions

""")
parser.add_argument("reference_genome", help="The reference genome to scan for the particular sequence")
parser.add_argument("sequence", help="The sequence from elsewhere to match")
args = parser.parse_args()
reference = args.reference_genome
print(reference[32:-16])
species_name = reference[32:-16]

database = reference[32:-16]+"/"+reference[32:-6]+'_db'
cmd = 'makeblastdb -in %s -dbtype nucl -out %s' % (reference, database)
subprocess.getoutput(cmd)

outfmt6 = reference[32:-16]+"/"+reference[32:-6]+'.outfmt6'
cmd2 = 'blastn -query %s -db %s -evalue=100000 -task "blastn" -outfmt 6 > %s' % (args.sequence, database, outfmt6)
subprocess.getoutput(cmd2)

extract = pd.read_csv(outfmt6, sep="\t", header=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

extract_df = extract.loc[extract['evalue']==min(extract['evalue'])]
extract_df = extract_df.reset_index()

extract_bed = pd.DataFrame(columns=['chrom', 'chromStart', 'chromEnd'])
for i in range(0, len(extract_df)):
    if extract_df['sstart'][i] <= extract_df['send'][i]:
        extract_bed = extract_bed.append({'chrom': extract_df['sseqid'][i],'chromStart': (extract_df['sstart'][i])-1, 'chromEnd': extract_df['send'][i]}, ignore_index=True)
    else:
        extract_bed = extract_bed.append({'chrom': extract_df['sseqid'][i],'chromStart': (extract_df['send'][i])-1, 'chromEnd': extract_df['sstart'][i]}, ignore_index=True)
extract_bed.sort_values(['chromStart'])


bedfile = reference[32:-16]+"/"+reference[32:-6]+'.ITS.bedfile'

extract_bed.to_csv(bedfile, sep='\t', header=False, index=False)

bedoutput = reference[32:-16]+"/"+reference[32:-6]+'.ITS.bedoutput.fasta'
cmd3 = 'bedtools getfasta -fo %s -fi %s -bed %s' % (bedoutput, reference, bedfile)
subprocess.getoutput(cmd3)

f = open(bedoutput, 'r')
filedata = f.read()
f.close()
new_data = filedata.replace('>', '>%s:' % species_name)
f = open(bedoutput, 'w')
f.write(new_data)
f.close()
