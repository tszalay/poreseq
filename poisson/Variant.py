import numpy as np
import pdb
import sys
from Bio import SeqIO

from LoadData import LoadAlignedEvents
from Util import RegionInfo



def Variant(ref_fasta, bamfile, fast5dir, var_fasta=None, muts=None, region=None, params={}, verbose=0):
    """Run variant-calling given required info.
    
    This function is the main one called by the setuptools entry points for
    variant calling, that is, 
        main() -> variant() -> Variant(...) -> PoissAlign.stuff(...)
        
    If you want to use variant calling from within Python, it's probably much
    easier just to do it by hand; copy the bottom part of this function.
    
    Args:
        ref_fasta (string): reference fasta file (alignment reference)
        bamfile (string): BAM-formatted file of alignments to reference
        fast5dir (string): folder where fast5 files are contained
        var_fasta (string): fasta sequences to compare to reference
        muts (list(MutInfo)): list of mutations to test, empty to test all points
        region (string): region string, or None for whole reference
        params (param dict): parameters to use for loading
        verbose (int): 0 for silent, 1 for steps, 2 for full mutations
        
    Returns:
        variant scores, if run with var_fasta, or MutScores if run
        with muts.
    """
    
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
    if muts is not None:
        if verbose > 0:
            sys.stderr.write("Variant calling {} using {} events\n".format(region,len(pa.events)))
            
        # offset the start to where pa.sequence begins
        for m in muts:
            m.start -= reginfo.start
            
        # now score mutations, either specified, or all of them
        if len(muts) > 0:
            mutscores = pa.ScoreMutations(muts)
        else:
            mutscores = pa.ScorePoints()
        for ms in mutscores:
            # now set the start back
            ms.start += reginfo.start
            sys.stdout.write(str(ms) + '\n')
            
        return mutscores
