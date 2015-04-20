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


def split_regions(fastafile, region_length, nfiles=None, perfile=None):

    refs = SeqIO.index(fastafile, "fasta")
    
    if nfiles is None and perfile is None:
        return

    fastabase = os.path.splitext(fastafile)[0]

    # now generate list of region strings
    regions = []
    for refid in refs:
        # load the reference sequence
        refseq = refs[refid]
        # check if we need to take sub-slices
        dl = region_length - 1000
        istart = 0
        iend = min(region_length,len(refseq))
        # step through max_length-1000 at a time (leave nice overlap)
        while istart < iend:
            regions.append('{}:{}:{}'.format(refid,istart,iend))
            iend = min(iend+dl,len(refseq))
            istart = min(istart+dl,len(refseq))

    # and write them to files
    if nfiles is not None:
        # split regions into a number of files
        chunks = []
        # open each file, write to them in parallel
        for i in range(nfiles):
            chunks.append(open(fastabase + '.{}.region'.format(i+1),'w'))
        # now loop through the regions
        for reg in regions:
            fileind = random.randint(0,nfiles-1)
            chunks[fileind].write(reg + '\n')
    elif perfile is not None:
        # or into files containing a number of regions
        fileind = -1
        f = None
        nwritten = perfile
        for reg in regions:
            # do we need to go on to the next file?
            if nwritten >= perfile:
                fileind += 1
                f = open(fastabase + '.{}.region'.format(fileind+1),'w')
                nwritten = 0
            f.write(reg + '\n')
            nwritten += 1
