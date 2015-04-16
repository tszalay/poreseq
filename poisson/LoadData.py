import os
import pysam
import pdb
import numpy as np
from Bio import SeqIO

from EventData import PoissEvent

def LoadReference(fastafile, refname=None):
    refs = SeqIO.index(fastafile, "fasta")
    refnames = list(refs.keys())
    if refname is None:
        if len(refs) == 1:
            refname = refnames[0]
        else:
            raise Exception("Multiple references in fasta, must specify one")
            
    return refs[refname].seq

def EventsFromBAM(eventdir, bamfile, reginfo, params):
    
    # first, load a bam file
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    
    # find the reference, if not specified
    if reginfo.name is None:
        if bamfile.nreferences > 1:
            raise Exception('Multiple references in BAM, one must be specified!')
        reginfo.name = bamfile.references[0]
    
    # and find which events correspond to a region
    bamevents = [x for x in bamfile.fetch(reference=reginfo.name,start=reginfo.start,end=reginfo.end)]
    
    # filter events by degree of overlap
    if 'min_overlap' in params:
        bamevents = [x for x in bamevents if x.get_overlap(reginfo.start,reginfo.end) >= params['min_overlap']]

    # and sort them by the amount of overlap, descending
    bamevents.sort(key = lambda x: x.get_overlap(reginfo.start,reginfo.end), reverse=True)
        
    # filter the event names to be unique
    #evnames = set([x.query_name for x in bamevents])
    
    # now filter to get one unique alignment per read
    #bamevents = [next(x for x in bamevents if x.query_name == evname) for evname in evnames]
    
    # and if a max number is set, take a subselection
    # but not randomly, use the most overlapping reads
    if 'min_coverage' in params and len(bamevents) < params['min_coverage']:
        raise Exception('Insufficient coverage!')
    if 'max_coverage' in params and len(bamevents) > params['max_coverage']:
        bamevents = bamevents[0:int(params['max_coverage'])]
    
    # now loop through and load
    events = []
    for bamev in bamevents:
        # the filename
        evfile = os.path.join(eventdir,bamev.query_name)
        
        # get the alignments
        ap = bamev.get_aligned_pairs()
        aps = np.array([x for x in ap if x[0] is not None and x[1] is not None])
        # offset the indices by the hard clip at the start, if present
        cig0 = bamev.cigar[0]
        if cig0[0] == 5:
            aps[:,0] += cig0[1]
        # and subtract off the start index from ref. side
        if reginfo.start > 0:
            aps[:,1] -= reginfo.start
            
        # get the reads, template and complement, flip if necessary
        for loc in ['t','c']:
            try:
                ev = PoissEvent(evfile,loc)
                if bamev.is_reverse:
                    ev.flip()
                ev.mapaligns(aps)
                events.append(ev)
            except Exception as e:
                print str(e.message)
                
    if not events:
        raise Exception('No aligned reads found!')
                
    return events
