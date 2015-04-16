#! /usr/bin/env python

import sys
import os
import random
from Bio import SeqIO

def split_fasta(fastafile, nchunks):

    refs = SeqIO.index(fastafile, "fasta")
    
    chunks = []
    
    fastabase = os.path.splitext(fastafile)[0]
    # open each file, write to them in parallel
    for i in range(nchunks):
        chunks.append(open(fastabase + '.{}.fasta'.format(i),'w'))
    
    # now loop through the sequences
    for ref in refs:
        fileind = random.randint(0,nchunks-1)
        chunks[fileind].write(refs[ref].format('fasta'))
    
def run():
    
    if len(sys.argv) < 3:
        print "usage: poisssplit fasta num_chunks"
        sys.exit()
            
    split_fasta(sys.argv[1],int(sys.argv[2]))
