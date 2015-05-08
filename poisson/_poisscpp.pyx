#distuils: language = c++

import cython
import numpy as np
cimport numpy as np
import random
import copy

from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string

from Util import MutationInfo,MutationScore

cdef extern from "cpp/AlignUtil.h":
    cdef cppclass AlignParams:
        AlignParams()
        int verbose 
        int realign_width
        int scoring_width
        double lik_offset

    cdef cppclass MutInfo:
        int start
        string orig
        string mut
    cdef cppclass MutScore:
        int start
        string orig
        string mut
        double score


cdef extern from "cpp/Sequence.h":
    cdef cppclass Sequence:
        string bases
        vector[int] states
        
        Sequence()
        Sequence(string seq)

cdef extern from "cpp/EventData.h":
    cdef cppclass ModelData:
        void setData(double* levmean, double* levstdv, double* sdmean, double* sdstdv, bool complement)
        void setParams(double prob_skip, double prob_stay, double prob_extend, double prob_insert)
    
    cdef cppclass EventData:
        ModelData model
        Sequence sequence
        vector[double] mean
        vector[double] ref_align
        vector[double] ref_like

        EventData()
        void setData(int length, double* mean, double* stdv, double* ref_align, double* ref_like)
        
cdef extern from "cpp/AlignData.h":
    cdef cppclass AlignData:
        Sequence sequence
        vector[EventData] events
        AlignParams params

cdef extern from "cpp/Mutations.h":
    cdef vector[MutInfo]    FindMutations(AlignData& data, vector[Sequence] seqs)
    cdef vector[MutInfo]    FindPointMutations(AlignData& data)

    cdef int                MakeMutations(AlignData& data, vector[MutScore] muts)
    cdef vector[MutScore]   ScoreMutations(AlignData& data, const vector[MutInfo]& muts)
    cdef vector[double]     ScoreAlignments(AlignData& data, double* likes)
        
cdef extern from "cpp/swlib.h":
    cdef cppclass SWAlignment:
        int score
        double accuracy
        vector[int] inds1
        vector[int] inds2
        
    cdef SWAlignment swfull(const string& seq1, const string& seq2)
    
cdef extern from "cpp/Viterbi.h":
    cdef vector[Sequence]   ViterbiMutate(vector[EventData]& events, int nkeep, 
                                double skip_prob, double stay_prob, 
                                double mut_min, double mut_max, bool verbose)
    
            
cdef double* getPr(np.ndarray[double, ndim=1, mode="c"] nparr):
    # get the pointer to a numpy ndarray, while ensuring it is contiguous
    return &nparr[0]

cdef vector[Sequence] PythonToSequences(pyseqs):
    # convert a list of strings to a C++ vector of Sequence classes
    cdef vector[Sequence] seqs
    
    for pyseq in pyseqs:
        seqs.push_back(Sequence(pyseq))
        
    return seqs

cdef vector[EventData] PythonToEvents(pyevents):
    # create C++ events from python PoissEvent objects
    cdef vector[EventData] events
    cdef EventData event
    
    cdef double* pr
    
    for pyev in pyevents:
        # make sure all the data is right
        pyev.makecontiguous()
        # create the event class, with default constructor
        event = EventData()
        # then set its data with above
        event.setData(pyev.mean.size,getPr(pyev.mean),getPr(pyev.stdv),getPr(pyev.ref_align),getPr(pyev.ref_like))
        # and its model's data
        event.model.setData(getPr(pyev.model.level_mean),getPr(pyev.model.level_stdv),
                            getPr(pyev.model.sd_mean),getPr(pyev.model.sd_stdv),pyev.model.complement)
        # also the model's params
        event.model.setParams(pyev.model.prob_skip,pyev.model.prob_stay,pyev.model.prob_extend,
                              pyev.model.prob_insert)
        # and set the sequence
        event.sequence = Sequence(pyev.sequence)
        
        n = 0
        for i in range(pyev.ref_align.size):
            if event.ref_align[i] > 0:
                n = event.ref_align[i]
        
        events.push_back(event)
        
    return events

cdef UpdatePythonEvents(pyevents, AlignData& data):
    '''Saves the computed values from from the C++ events to Python events
    (ref_align and ref_likes)'''
    for i,ev in enumerate(pyevents):
        ev.ref_align[:] = data.events[i].ref_align
        ev.ref_like[:] = data.events[i].ref_like
    return pyevents
    
cdef AlignData PythonToAlignData(obj):
    '''Creates C++ AlignData from Python's PoissAlign'''
    cdef AlignData data
    data.sequence = Sequence(obj.sequence)
    data.events = PythonToEvents(obj.events)
    if 'verbose' in obj.params:
        data.params.verbose = obj.params['verbose']
    if 'lik_offset' in obj.params:
        data.params.lik_offset = obj.params['lik_offset']
    if 'realign_width' in obj.params:
        data.params.realign_width = obj.params['realign_width']
    if 'scoring_width' in obj.params:
        data.params.scoring_width = obj.params['scoring_width']
            
    return data
    
def swalign(seq1,seq2):
    """Smith-Waterman align two sequences.
    
    Uses a full-matrix SW alignment function to calculate alignments of seq1
    and seq2. This is a wrapper for C++, so it is fast, but expect the memory
    to go as N^2; in practice, sequences longer than 20kb become troublesome.
    
    Args:
        seq1 (string): DNA sequence to align
        seq2 (string): DNA sequence to align
    
    Returns:
        Tuple of (accuracy, pairs): returns the accuracy (in %) and a list of
        tuples of pairwise aligned indices.
    """
    cdef SWAlignment align = swfull(seq1,seq2)
    pairs = []
    for i in range(align.inds1.size()):
        pairs.append((align.inds1[i],align.inds2[i]))
    return (align.accuracy, pairs)

class PoissAlign:
    """Class containing all data associated with reads aligned to a reference.

    This is the go-to class for loading and analyzing aligned event data.
    All class methods are in-place and mutable/non-const, use .Copy() if non-
    destructive behavior is desired.
    
    Note also that loaded events aren't realigned by default, until any of the
    .Score... functions are called.

    Example:

        params = poisson.LoadParams('./params.conf')
        reginfo = poisson.RegionInfo('10000:20000')
        # pa is the PoissAlign class
        pa = poisson.LoadAlignedEvents('reference.fasta','alignment.bam',
                                       '/media/run-data',reginfo,params)
        scores = pa.ScoreEvents()  # this forces realignment
        print pa.events[2].ref_align[2000:3000]    
        
    Attributes:
        sequence (string): the reference that events are currently aligned to
        events (list(PoissEvent)): aligned PoissEvent classes
        params (params dict): parameters to use for sub-functions
    """
    
    def __init__(self):
        '''Empty initializer'''
        self.sequence = ""
        self.events = []
        self.params = {}
        
    def Copy(self):
        '''Wrapper around copy.deepcopy'''
        return copy.deepcopy(self)
        
    def Coverage(self):
        """Calculate depth of coverage across the reference.
    
        Returns:
            double x N ndarray: the number of events aligned at each of the
                N bases of self.sequence
        """
        cov = np.zeros(len(self.sequence))
        for ev in self.events:
            nzs = ev.ref_align[ev.ref_align>0]
            minind = nzs[0]
            maxind = np.minimum(nzs[-1],len(cov)-1)
            cov[minind:maxind] += 1
            
        return cov
        
    def RealignTo(self, newseq):
        """Realign all of the events to a new reference sequence using swalign.
    
        This function uses Smith-Waterman to in-place realign the loaded events to a
        new reference (and updates self.sequence). Useful for scoring a similar
        but not identical sequence.
        
        Args:
            newseq (string): sequence that aligns to self.sequence
        
        Returns:
            None
        """
        # calculate the alignments
        align = swalign(self.sequence, newseq)
        if align[0] < 0.6:
            raise Exception('Error rate too large for realignment!')
        # now actually map the alignments
        [x.mapaligns(np.array(align[1])) for x in self.events]
        # and set the sequence
        self.sequence = newseq
        
    def ScoreEvents(self):
        """Calculate likelihood score of all aligned events.
    
        This function realigns all events to the reference and recalculates
        all of the total likelihood scores for each event.
        
        Returns:
            M x double of likelihood scores for each event in self.events
        """
        # calculate likelihoods over all of the strands
        cdef AlignData data = PythonToAlignData(self)
        # and get the scores, without the aligned likelihoods
        cdef vector[double] scores = ScoreAlignments(data,NULL)
        return scores
        
    def ScorePoints(self):
        """Calculate mutation score of all single base mutations.
    
        This function returns a list of MutationScore() objects for each
        possible single base mutation in the sequence. Trivial mutations
        (eg. A -> A) are omitted, but redundant ones (eg. deletions in
        homopolymer regions) are not.
        
        Returns:
            list(MutationScore) for all single-base mutations
        """
        
        # get the score of all single-base mutations
        cdef AlignData data = PythonToAlignData(self)
        # override scoring width with point mutation width if available
        if 'point_width' in self.params:
            data.params.scoring_width = self.params['point_width']
            
        cdef vector[MutInfo] mutations = FindPointMutations(data)
        cdef vector[MutScore] scores = ScoreMutations(data,mutations)
        
        pyscores = []
        for i in range(scores.size()):
            s = MutationScore()
            s.start = scores[i].start
            s.orig = scores[i].orig
            s.mut = scores[i].mut
            s.score = scores[i].score
            pyscores.append(s)
        
        return pyscores
        
    def ScoreMutations(self, muts):
        """Calculate mutation score of specific mutations.
    
        This function returns a list of MutationScore() objects for each
        MutationInfo in muts, with the score calculated.
        
        Args:
            muts (list(MutationInfo)): list of mutations to test
        
        Returns:
            list(MutationScore) for all mutations given
        """
        # get the score of all single-base mutations
        cdef AlignData data = PythonToAlignData(self)
            
        cdef vector[MutInfo] mutations
        cdef MutInfo mi
        
        for m in muts:
            mi.start = m.start
            mi.orig = m.orig
            mi.mut = m.mut
            mutations.push_back(mi)
            
        cdef vector[MutScore] scores = ScoreMutations(data,mutations)
        
        pyscores = []
        for i in range(scores.size()):
            s = MutationScore()
            s.start = scores[i].start
            s.orig = scores[i].orig
            s.mut = scores[i].mut
            s.score = scores[i].score
            pyscores.append(s)
        
        return pyscores
        
    def ApplyMuts(self, pymuts):
        """Make specific mutations.
    
        Don't use this function, its behavior is a bit dodgy. The C++ code 
        will attempt to make all positive-scoring mutations, but if there are
        ones too close together, they will be saved for a second pass and then
        re-scored and re-tested. This function should really just make them
        forcibly.
        
        Args:
            muts (list(MutationInfo)): list of mutations to make        
        """
        # get the score of all single-base mutations
        cdef AlignData data = PythonToAlignData(self)
        if 'point_width' in self.params:
            data.params.scoring_width = self.params['point_width']
            
        cdef vector[MutScore] scores
        cdef MutScore ms
        for m in pymuts:
            ms.start = m.start
            ms.orig = m.orig
            ms.mut = m.mut
            ms.score = m.score
            scores.push_back(ms)
            
        MakeMutations(data,scores)
        self.sequence = data.sequence.bases
        self.events = UpdatePythonEvents(self.events, data)
        
        
    def Mutate(self,seqs='self',reps=4):
        """Take similar sequences and use them to mutate the consensus sequence.
    
        This function is the first of two used in consensus error correction.
        It can either take all of the 2D sequences in loaded events, a user-
        supplied list of sequences, or generate Viterbi sequences, and it uses
        these sequences as seeds for mutations to test and make to improve the
        overall accuracy.
        
        The mutation finding is cached in AlignData for later iterations, so
        even though it 'finds' mutations multiple times, only the first one
        takes a long time.
        
        Example, from overlap alignment:
        
            params = poisson.LoadParams('./params.conf')
            reginfo = poisson.RegionInfo('ch_22_strand_55.fast5:10000:20000')
            # pa is the PoissAlign class
            pa = poisson.LoadAlignedEvents('extracted.fasta','alignment.bam',
                                       '/media/run-data',reginfo,params)
            pa.Mutate()
       
        Args:
            seqs (string, list(string)): either 'self','viterbi' or a list of
                            sequences to use
            reps (int): how many times to repeat the mutation testing
        
        Returns:
            nbases (int): total number of bases mutated
        """
        cdef AlignData data = PythonToAlignData(self)
        
        cdef vector[Sequence] sequences
        
        if seqs == 'self':
            # use every other sequence to account for template/complement
            seqs = [x.sequence for x in self.events[::2]]
        elif seqs == 'viterbi':
            seqs = None
            sequences = ViterbiMutate(data.events,16,0.05,0.01,0.33,0.75,self.params['verbose'])
            
        if seqs:
            sequences = PythonToSequences(seqs)
        
        cdef vector[MutInfo] mutations
        cdef vector[MutScore] scores
        totbases = 0
        for i in range(reps):
            mutations = FindMutations(data,sequences)
            scores = ScoreMutations(data,mutations)
            nbases = MakeMutations(data,scores)
            if nbases == 0:
                break
            totbases += nbases
            
        self.sequence = data.sequence.bases
        self.events = UpdatePythonEvents(self.events, data)
        return totbases
        
    def Refine(self):
        """Test all single-base mutations to improve the consensus sequence.
    
        This function is the second of two used in consensus error correction.
        It brute force tests each single base substitution, deletion, and insertion
        in order to find and make any that improve the overall score.
        
        This is sped up by the use of the narrower point_width as specified in
        params.
        
        For example, to do basic variant calling:
        
            params = poisson.LoadParams('./params.conf')
            reginfo = poisson.RegionInfo('10000:20000')
            # pa is the PoissAlign class
            pa = poisson.LoadAlignedEvents('reference.fasta','alignment.bam',
                                       '/media/run-data',reginfo,params)
            pa.Refine()
            
        at which point pa.sequence is the new sequence with all variants,
        and can be compared to the original reference to find differences.
       
        Returns:
            nbases (int): total number of bases mutated
        """
        
        cdef AlignData data = PythonToAlignData(self)
        # override scoring width with point mutation width if available
        if 'point_width' in self.params:
            data.params.scoring_width = self.params['point_width']
        cdef vector[MutInfo] mutations = FindPointMutations(data)
        cdef vector[MutScore] scores = ScoreMutations(data,mutations)
        nbases = MakeMutations(data,scores)
        self.sequence = data.sequence.bases
        self.events = UpdatePythonEvents(self.events, data)
        return nbases

    