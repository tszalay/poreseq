#! /usr/bin/env python

import sys
import os
import random
from Bio import SeqIO

def split_fasta(fastafile, nchunks=None, nseqs=None):

    refs = SeqIO.index(fastafile, "fasta")
    
    if nchunks is None and nseqs is None:
        return

    fastabase = os.path.splitext(fastafile)[0]

    if nchunks is not None:
        # split seqs into a number of files
        chunks = []
        # open each file, write to them in parallel
        for i in range(nchunks):
            chunks.append(open(fastabase + '.{}.fasta'.format(i+1),'w'))
        # now loop through the sequences
        for ref in refs:
            fileind = random.randint(0,nchunks-1)
            chunks[fileind].write(refs[ref].format('fasta'))
    elif nseqs is not None:
        # or into files containing a number of sequences
        fileind = -1
        f = None
        nwritten = nseqs
        for ref in refs:
            # do we need to go on to the next file?
            if nwritten >= nseqs:
                fileind += 1
                f = open(fastabase + '.{}.fasta'.format(fileind+1),'w')
                nwritten = 0
            f.write(refs[ref].format('fasta'))
            nwritten += 1
