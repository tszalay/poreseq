import numpy as np
import pdb
import sys
from Bio import SeqIO

from LoadData import LoadAlignedEvents
from Util import RegionInfo



def Variant(ref_fasta, bamfile, fast5dir, var_fasta, region=None, params={}, verbose=0):

    pa = LoadAlignedEvents(ref_fasta,bamfile,fast5dir,RegionInfo(region),params)
    pa.params['verbose'] = verbose
    
    # load variant fasta files
    variants = SeqIO.index(var_fasta, "fasta")

    if verbose > 0:
        sys.stderr.write("Variant calling {} variants with {} bases using {} events\n".format(len(variants),len(pa.sequence),len(pa.events)))

    # go through each variant, re-align and call
    variantscores = {}
    for variant in variants:
        # create a copy so we can realign it
        pav = pa.Copy()
        varseq = str(variant.seq)
        pav.RealignTo(varseq)
        # ok, now it's aligned to the variant, get the score over all strands
        score = np.sum(pav.ScoreEvents())
        print score

    return pa.sequence
