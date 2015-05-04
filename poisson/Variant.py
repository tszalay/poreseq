import numpy as np
import pdb
import sys
from Bio import SeqIO

from LoadData import LoadAlignedEvents
from Util import RegionInfo



def Variant(ref_fasta, bamfile, fast5dir, var_fasta=None, muts=[], region=None, params={}, verbose=0):

    # load all of the aligned events, same as for consensus
    reginfo = RegionInfo(region)
    pa = LoadAlignedEvents(ref_fasta,bamfile,fast5dir,reginfo,params)
    pa.params['verbose'] = verbose
    
    # variant call if sequences specified directly
    if var_fasta is not None:
        variants = SeqIO.index(var_fasta, "fasta")
    
        if verbose > 0:
            sys.stderr.write("Variant calling {} variant sequences with {} bases using {} events\n".format(len(variants),len(pa.sequence),len(pa.events)))
    
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
            variantscores[record.id] = dscore
            
        return variantscores
    
    # otherwise, score the mutations as specified
    if len(muts) > 0:        
        if verbose > 0:
            sys.stderr.write("Variant calling {} using {} events\n".format(region,len(pa.events)))

        # offset the start to where pa.sequence begins
        for m in muts:
            m.start -= reginfo.start
        # now check to see if they're within 
        # and score them
        mutscores = pa.ScoreMutations(muts)
        for ms in mutscores:
            # now set the start back
            ms.start += reginfo.start
            sys.stdout.write(str(ms) + '\n')
