import sys
import os
import random
from Bio import SeqIO

def split_fasta(fastafile, nchunks=None, nseqs=None):
    """Split all of the sequences in fastafile into multiple fasta files.
    
    This is one way to split up the error correction of many sequences (eg.
    fast5-extracted sequences) into sub-tasks. Specifying nchunks splits up 
    sequences into nchunks files, while nseqs puts nseqs sequences into each file.
    
    Args:
        fastafile (string): fasta file containing many sequences
        nchunks(int): how many files to split to
        nseqs(int): how many sequences to write per file    
    """
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


def split_regions(fastafile, region_length, nfiles=None, perfile=None, userefs=None):
    """Generates region strings for splitting up sequences in fastafile.
    
    This is the main function for splitting up work into overlapping chunks of
    size region_length. If nfiles or perfile is specified, it outputs the region
    strings to files, otherwise just returns them. For example:
    
    >sequence1
    ACCCGT....... (length 35000)
    
    will lead to the following region strings (for region_length=10000):
    
    sequence1:0:10000
    sequence1:9000:19000
    sequence1:18000:28000
    sequence1:27000:35000
    
    which are either returned or saved to a file. If userefs is not specified,
    this is done for all sequences in fastafile.
    
    Args:
        fastafile (string): reference fasta file
        region_length (int): length of each region
        nfiles(int): write region strings to N files
        perfile(int): write M region strings to each file
        userefs(list(string)): get regions only from these references
    
    Returns:
        list(string): list of region strings, if no file output specified
        None, if file output specified
    """
    
    refs = SeqIO.index(fastafile, "fasta")
    
    region_length = int(region_length)
    # now generate list of region strings
    regions = []
    for refid in refs:
        # if sub-refs to use specified and not there, skip
        if userefs is not None and refid not in userefs:
            continue
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
    
    # if we didn't specify how to split, we're just returning all of the
    # split regions as-is as a utility function
    if nfiles is None and perfile is None:
        return regions

    fastabase = os.path.splitext(fastafile)[0]

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
