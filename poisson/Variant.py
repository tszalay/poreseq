import numpy as np
import pdb
import sys
from Bio import SeqIO

from LoadData import LoadAlignedEvents
from Util import RegionInfo



def Variant(ref_fasta, bamfile, fast5dir, var_fasta, region=None, params={}, verbose=0):

    # load all of the aligned events, same as for mutation
    pa = LoadAlignedEvents(ref_fasta,bamfile,fast5dir,RegionInfo(region),params)
    pa.params['verbose'] = verbose
    
    # load variant fasta files
    variants = SeqIO.index(var_fasta, "fasta")

    if verbose > 0:
        sys.stderr.write("Variant calling {} variants with {} bases using {} events\n".format(len(variants),len(pa.sequence),len(pa.events)))

    basescore = np.sum(pa.ScoreEvents())

    # go through each variant, re-align and call
    variantscores = {}
    for variant in variants:
        # create a copy so we can realign it
        pav = pa.Copy()
        record = variants[variant]
        varseq = str(record.seq)
        pav.RealignTo(varseq)
        # ok, now it's aligned to the variant, get the score over all strands
        dscore = np.sum(pav.ScoreEvents()) - basescore
        sys.stdout.write('{}, {}\n'.format(record.id, dscore))
