import numpy as np
from LoadData import LoadAlignedEvents
from Util import RegionInfo
import poisscpp
import pdb
import sys

def Mutate(fastafile, bamfile, fast5dir, region=None, params={}, verbose=0, test=False, reps=4):
    """Run consensus-calling mutations given required info.
    
    This function is the main one called by the setuptools entry points for
    consensus calling, that is, 
        main() -> consensus() -> Mutate(...) -> PoissAlign.Mutate(...)
        
    Use this function if you want to basically run poisson consensus from
    within Python.
    
    It takes care of loading events and mutating in different ways, and if the test
    flag is specified, it starts with a low-accuracy 2D sequence instead of the ref.
    
    Args:
        fastafile (string): reference fasta file (alignment reference)
        bamfile (string): BAM-formatted file of alignments to reference
        fast5dir (string): folder where fast5 files are contained
        region (string): region string, or None for whole reference
        params (param dict): parameters to use for loading
        verbose (int): 0 for silent, 1 for steps, 2 for full mutations
        test (boolean): start with low-accuracy 2D sequence?
        reps (int): number of iterations for mutation and refinement
    
    Returns:
        tuple(sequence (string), acc (float)): mutated higher-accuracy consensus sequence
        along with its accuracy relative to the reference
    """

    if 'verbose' not in params:
        params['verbose'] = 0
    
    pa = LoadAlignedEvents(fastafile,bamfile,fast5dir,RegionInfo(region),params)
    
    # and the loaded reference sequence
    refseq = pa.sequence
    
    # test automatically sets verbose
    if test and verbose == 0:
        verbose = 1
    
    # we know our algorithm doesn't do great for 1 or 2 events
    # in which case we can just shortcut and return the starting seq
    if len(pa.events) < 5:
        if verbose > 0:
            sys.stderr.write("Coverage is 1 or 2, not mutating...\n")
        return (refseq, 100)

    if verbose > 0:
        sys.stderr.write("Mutating {} bases using {} events\n".format(len(refseq),len(pa.events)))

    # if test mode, pick a sequence from event sequences
    if test:
        seq = ""
        for ev in pa.events:
            pairs = poisscpp.swalign(ev.sequence,refseq)[1]
            if pairs[-1][1]-pairs[0][1] > len(seq):
                seq = ev.sequence[pairs[0][0]:pairs[-1][0]]
        pa.sequence = seq

    if test:
        sys.stderr.write("Starting accuracy: " + str(round(poisscpp.swalign(pa.sequence,refseq)[0],1)) + "%\n")

    pa.Mutate(reps=reps)
    
    if verbose>0:
        acc = poisscpp.swalign(pa.sequence,refseq)[0]
        sys.stderr.write("Accuracy: " + str(round(acc,1)) + "%\n")

    for i in range(reps):
        
        pa.Mutate(seqs='viterbi')
        nbases = pa.Refine()
        
        if verbose>0:
            acc = poisscpp.swalign(pa.sequence,refseq)[0]
            sys.stderr.write("Accuracy: " + str(round(acc,1)) + "%\n")
        if nbases == 0:
            break
        
    # trim ends as requested
    if 'end_trim' in params and len(pa.sequence) > 2*params['end_trim']:
        pa.sequence = pa.sequence[int(params['end_trim']):-int(params['end_trim'])]

        
    # find the aligned sequence stats
    acc,inds = poisscpp.swalign(pa.sequence,refseq)

    if verbose>0:
        errs = np.sum(np.array(inds)==0,0)
        sys.stderr.write("Final accuracy: " + str(round(acc,1)) + "%\n")
        sys.stderr.write("Insertions: {}, Deletions: {}\n".format(errs[0],errs[1]))
        sys.stderr.write("Final coverage: " + str(round(np.mean(pa.Coverage()),1)) + "X\n")

    return (pa.sequence, acc)
