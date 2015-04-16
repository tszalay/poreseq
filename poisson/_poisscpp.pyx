#distuils: language = c++

import cython
import numpy as np
cimport numpy as np
import random
import copy

from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string

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
    return &nparr[0]

cdef vector[Sequence] PythonToSequences(pyseqs):

    cdef vector[Sequence] seqs
    
    for pyseq in pyseqs:
        seqs.push_back(Sequence(pyseq))
        
    return seqs

cdef vector[EventData] PythonToEvents(pyevents):
    
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
    for i,ev in enumerate(pyevents):
        ev.ref_align[:] = data.events[i].ref_align
    return pyevents
    
cdef AlignData PythonToAlignData(obj):
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
    cdef SWAlignment align = swfull(seq1,seq2)
    pairs = []
    for i in range(align.inds1.size()):
        pairs.append((align.inds1[i],align.inds2[i]))
    return (align.accuracy, pairs)


class PoissAlign:
    def __init__(self):
        self.sequence = ""
        self.events = []
        self.params = {}
        
    def Copy(self):
        return copy.deepcopy(self)
        
    def Coverage(self):
        cov = np.zeros(len(self.sequence))
        for ev in self.events:
            nzs = ev.ref_align[ev.ref_align>0]
            minind = nzs[0]
            maxind = np.minimum(nzs[-1],len(cov)-1)
            cov[minind:maxind] += 1
            
        return cov

        
    def Mutate(self,seqs='self',reps=4):
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

    