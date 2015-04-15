import numpy as np
from LoadData import *
from ParamData import *
from Util import *
from poisson import poisscpp
from poisscpp import PoissAlign
import pdb
import sys

def DoMutation(fastafile, bamfile, fast5dir, region=None, paramfile=None, verbose=0, test=False):

    reginfo = RegionInfo(region)
    refseq = str(LoadReference(fastafile,reginfo.name))
    
    # set region indices to full fasta, if none specified
    if reginfo.start is None and reginfo.end is None:
        reginfo.start = 0
        reginfo.end = len(refseq)-1
        
    params = {}
    # load params, if specified
    if paramfile is not None:
        params = LoadParams(paramfile)
    
    # now load the events
    events = EventsFromBAM(fast5dir,bamfile,reginfo,params)

    if len(params) > 0:
        [x.setparams(params) for x in events]
    
    if verbose > 0:
        sys.stderr.write("Running with " + str(len(events)) + " events\n")

    #pdb.set_trace()
    refseq = refseq[reginfo.start:reginfo.end]

    if test:
        seq = ""
        for ev in events:
            pairs = poisscpp.swalign(ev.sequence,refseq)[1]
            if pairs[-1][1]-pairs[0][1] > len(seq):
                seq = ev.sequence[pairs[0][0]:pairs[-1][0]]
    else:
        seq = refseq
        
    pa = PoissAlign()
    pa.sequence = seq
    pa.events = events
    pa.params = params
    pa.params['verbose'] = verbose

    if test:
        sys.stderr.write("Starting accuracy: " + str(round(poisscpp.swalign(pa.sequence,refseq)[0],1)) + "%\n")

    pa.Mutate()
    
    if test:
        acc = poisscpp.swalign(pa.sequence,refseq)[0]
        sys.stderr.write("Accuracy: " + str(round(acc,1)) + "%\n")

    for i in range(3):
        pa.Mutate(seqs='viterbi')
        nbases = pa.Refine()
        
        if test:
            acc = poisscpp.swalign(pa.sequence,refseq)[0]
            sys.stderr.write("Accuracy: " + str(round(acc,1)) + "%\n")
        if nbases == 0:
            break
        
    # trim ends as requested
    if 'end_trim' in params and len(pa.sequence) > 2*params['end_trim']:
        pa.sequence = pa.sequence[params['end_trim']:-params['end_trim']]

        
    # find the aligned sequence stats
    if test:
        acc,inds = poisscpp.swalign(pa.sequence,refseq)
        errs = np.sum(np.array(inds)==0,0)
        sys.stderr.write("Insertions: {}, Deletions: {}\n".format(errs[0],errs[1]))
        return (pa.sequence, acc)
    
    return pa.sequence