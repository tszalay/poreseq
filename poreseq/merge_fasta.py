import sys
import os
from Bio import SeqIO
from Util import RegionInfo
from poreseqcpp import swalign
import pdb

def merge_seqs(seq1, seq2, overlap):
    '''Merge seq1 and seq2 assuming there are 'overlap' bases in common.

    Aligns the last 'overlap' bases of seq1 with the first 'overlap' bases of
    seq2, and if they align well enough, it joins them halfway.

    Args:
        seq1 (string): first sequence, goes on the left
        seq2 (string): second sequence, on the right
        overlap (int): usually 1000
        
    Returns:
        Joined sequence
    '''
    # assumes the last "overlap" bases of seq1 and the first "overlap" bases of seq2 overlap
    i0 = -overlap
    i1 = overlap
    if len(seq1) < overlap:
        i0 = 0
    if len(seq2) < overlap:
        i1 = len(seq2)-1

    acc,inds = swalign(seq1[i0:],seq2[:i1])
    # put some constraints on the accuracy
    if acc < 0.70:
        raise Exception('Insufficient accuracy for overlap')
    # now put the break point halfway between
    inds = [x for x in inds if x[0]>0 and x[1]>0]
    imid = inds[int(len(inds)/2)]
    i0 += imid[0]
    i1 = imid[1]
    return seq1[:i0] + seq2[i1:]
    

def merge_fasta(fastafiles, fastaout):
    """Merge all of the sequences in fastafiles and put them in fastaout.
    
    Loads all of the sequences in fastafiles. If any of the headers have
    positional regions appended (eg. >filename.fast5:10000:20000), uses
    these regions to attempt to merge all of the sequences and save the full
    >filename.fast5 as a single sequence.
    
    Args:
        fastafiles (list(string): list of fasta filenames
        fastaout (string): output fasta filename
    """
    fragments = {}

    # loop through each given fasta file and load them
    for fasta in fastafiles:
        refs = SeqIO.index(fasta, "fasta")
        # go through each ref header in this file
        for ref in refs:
            # this header is a combination of the filename and the header
            # so let's parse it
            reg = RegionInfo(ref)
            # and put all corresponding pieces in a list, along with the sequence
            if reg.name not in fragments:
                fragments[reg.name] = []
            fragments[reg.name].append((reg,str(refs[ref].seq)))

    # now fragments contains hopefully all the pieces, let's put them together
    # one sequence at a time

    outfile = open(fastaout,'w')

    for ref in fragments:
        seqlist = fragments[ref]
        # first, sort the list in order of start index
        seqlist.sort(key=lambda x: x[0].start)
        # and now reduce down to a single sequence
        seq = reduce(lambda x,y: merge_seqs(x,y,1000), [x[1] for x in seqlist])
        outfile.write('>{}\n{}\n'.format(ref,seq))

    
