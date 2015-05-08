import os
import pysam
import pdb
import numpy as np
from Bio import SeqIO

from EventData import PoissEvent
from poisson import poisscpp

def LoadAlignedEvents(fastafile, bamfile, eventdir, reginfo, params):
    """Load all of the events aligned to a reference using a bam file.
    
    This is the main function for loading events subject to various parameters
    and constraints. The function returns a PoissAlign object defined in the
    pyx file, which can then be used for subsequent processing.
    
    The params argument is particularly important here, for min_coverage,
    max_coverage, min_overlap, etc.
    
    Args:
        fastafile (string): reference fasta file (alignment reference)
        bamfile (string): BAM-formatted file of alignments to reference
        eventdir (string): folder where fast5 files are contained
        reginfo (RegionInfo): which region to load? RegionInfo(None) for all
        params (param dict): parameters to use for loading
    
    Returns:
        PoissAlign object containing reference sequence, aligned events, and params
    """
    # load the reference sequence    
    refseq = str(LoadReference(fastafile,reginfo.name))    
    # set region indices to full fasta, if none specified
    if reginfo.start is None and reginfo.end is None:
        reginfo.start = 0
        reginfo.end = len(refseq)       
    # now load the events
    events = EventsFromBAM(eventdir,bamfile,reginfo,params)

    # set the parameters for the loaded events
    if len(params) > 0:
        [x.setparams(params) for x in events]
    
    # take specified region of refseq
    refseq = refseq[reginfo.start:reginfo.end]
    
    # create poissalign object with loaded settings
    pa = poisscpp.PoissAlign()
    pa.sequence = refseq
    pa.events = events
    pa.params = params
    
    return pa

def LoadReference(fastafile, refname=None):
    '''Loads a reference sequence from a fasta file'''
    
    refs = SeqIO.index(fastafile, "fasta")
    refnames = list(refs.keys())
    if refname is None:
        if len(refs) == 1:
            refname = refnames[0]
        else:
            raise Exception("Multiple references in fasta, must specify one")
            
    return refs[refname].seq

def EventsFromBAM(eventdir, bamfile, reginfo, params):
    """Does the heavy lifting of LoadAlignedEvents, returns list of events.
    
    Main steps:
    
        - Find events aligned to specified region in the BAM file
        - Sort these events by amount of overlap
        - Keep unique events up to max_coverage
        - Load events (temp. and comp.) from fast5 files
        - Map alignments from events' own 2D seqs to the reference sequence
        - Return a list of all loaded and aligned events
    """
    
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
    
    # now go through and collect the events we want to keep (up to max coverage, and unique)
    bamnames = []
    newevents = []
    for bamev in bamevents:
        if not bamev.query_name in bamnames:
            bamnames.append(bamev.query_name)
            newevents.append(bamev)
        # stop if we have hit our desired coverage level
        if 'max_coverage' in params and len(newevents) >= params['max_coverage']:
            break
        
    bamevents = newevents
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
