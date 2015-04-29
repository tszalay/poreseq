import h5py
import pdb
import numpy as np
import copy
from Bio.Seq import Seq

def MakeContiguous(pyobj):
    attrs = [attr for attr in dir(pyobj) if not callable(attr) and not attr.startswith("__")]
    for attr in attrs:
        val = getattr(pyobj,attr)
        if type(val) == np.ndarray:
            setattr(pyobj,attr,np.ascontiguousarray(val,dtype='f8'))

class PoissModel:
    def __init__(self):
        self.level_mean = []
        self.level_stdv = []
        self.sd_mean = []
        self.sd_stdv = []
        self.prob_skip = 0.1
        self.prob_stay = 0.1
        self.prob_extend = self.prob_stay
        self.prob_insert = 0.01
        self.name = ''
        self.complement = False

class PoissEvent():
    def __init__(self, filename,typ):
        # load fast5 file
        f = h5py.File(filename,'r')
        # and try to read the data from the location
        loc = 'template'
        if typ[0] == 'c':
            loc = 'complement'
        evdata = f['/Analyses/Basecall_2D_000/BaseCalled_'+loc+'/Events']
        modeldata = f['/Analyses/Basecall_2D_000/BaseCalled_'+loc+'/Model']
        attdata = f['/Analyses/Basecall_2D_000/Summary/basecall_1d_'+loc].attrs
        
        # also get the Oxford-called 2D sequence
        seqdata = f['/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq']
        # load sequence (parse from Fastq)
        self.sequence = seqdata[()].split('\n')[1]

        # and the corresponding alignment to the reference
        aldata = f['/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment']
        
        # now initialize the alignment
        alinds = aldata[loc]
        kmers = aldata['kmer']
        seqinds = 0*alinds
        curind = 0
        for i in range(len(aldata[loc])):
            curind = self.sequence.find(kmers[i],curind)
            seqinds[i] = curind

        shift = attdata['shift']
        scale = attdata['scale']
        scalesd = attdata['scale_sd']
        drift = attdata['drift']
        var = attdata['var']
        varsd = attdata['var_sd']
        
        # load and scale event data
        self.mean = evdata['mean']
        self.stdv = evdata['stdv']
        self.length = evdata['length']
        self.start = evdata['start']
        self.mean = self.mean - drift*(self.start - self.start[0])
        self.ref_align = 0*self.mean
        self.ref_like = 0*self.stdv
        
        # now seed ref_aligns to the self-alignment
        # first, take everything from 'alignment' that is aligned
        lvlinds = alinds > 0
        # now set those indices in ref_align
        self.ref_align[alinds[lvlinds]] = seqinds[lvlinds]
        
        # load and scale model data
        model = PoissModel();
        model.level_mean = modeldata['level_mean']*scale + shift
        model.level_stdv = modeldata['level_stdv']*var
        model.sd_mean = modeldata['sd_mean']*scalesd
        model.sd_stdv = modeldata['sd_stdv']/np.sqrt(varsd)
        model.name = attdata['model_file']
        model.complement = (loc=='complement')
        self.model = model
        
        # now if it's a complement event, we need to flip it, but leave seq
        if self.model.complement:
            self.flip(False)
    
    def copy(self):
        return copy.deepcopy(self)        
        
    def flip(self, flip_sequence=True):
        
        # flip event members
        self.mean = self.mean[::-1]
        self.stdv = self.stdv[::-1]
        self.length = self.length[::-1]
        self.start = self.start[::-1]
        self.ref_align = self.ref_align[::-1]
        self.ref_like = self.ref_like[::-1]
        
        # and model members, by calculating flipped bit indices
        # 1023-x takes care of complementing
        flips = 1023-np.arange(1024)
        # now we just need to reverse
        flips = (((flips&0b11)<<8) | ((flips>>8)&0b11) | ((flips&0b1100)<<4) | 
                    ((flips>>4)&0b1100) | (flips&0b110000))
        
        # then just flip model props
        self.model.level_mean = self.model.level_mean[flips];
        self.model.level_stdv = self.model.level_stdv[flips];
        self.model.sd_mean = self.model.sd_mean[flips];
        self.model.sd_stdv = self.model.sd_stdv[flips];
        
        if flip_sequence:
            # now flip sequence if necessary, keeping it a string
            self.sequence = str(Seq(self.sequence).reverse_complement())
            # and shuffle the ref_align indices too, but only the aligned ones
            ra0 = self.ref_align > 0
            self.ref_align[ra0] = len(self.sequence) - self.ref_align[ra0]
            
        self.makecontiguous()
            
    def mapaligns(self, pairs):
        # re-map ref_aligns using provided pairs
        # first index is to current (ref_align) sequence, second is to new ref
        # (should be numpy array)

        # clear ref_align
        refal = self.ref_align
        ra0 = refal > 0
        self.ref_align = 0*self.ref_align
        # now make the pairs unique in x (own seq)
        (_,uinds) = np.unique(pairs[:,0],return_index=True)
        pairs = pairs[uinds,:]
        # and now use interpolation to remap, return -1 out of range
        self.ref_align[ra0] = np.round(np.interp(refal[ra0],pairs[:,0],pairs[:,1],0,0))
        
        # make sure all of our arrays are in order
        self.makecontiguous()
        
    def makecontiguous(self):
        MakeContiguous(self)
        MakeContiguous(self.model)
    
    def getrefstats(self):
        # figure out statistics of alignment to ref_align (skips,stays,etc)
        bins = np.bincount(np.int64(self.ref_align[self.ref_align>=0]))
        # ref_align positions not appearing
        skips = np.sum(bins[1:] == 0)
        # ref_align positions appearing more than once, count total
        stays = np.sum(np.maximum(0,bins[1:]-1))
        # and insertions
        inserts = np.sum(self.ref_align < 0)
        # and total number aligned
        total = float(np.sum(self.ref_align != 0))
        return (skips/total,stays/total,inserts/total)
        
    def setparams(self, params):
        # set skip/stay params
        # expects keys of the form 'insert_t' and sets properties
        # of the form 'prob_insert' in eg. template models only
        for k in params:
            paramname = 'prob_' + k[:-2]
            
            if paramname not in dir(self.model):
                continue
            if ((k[-2:] == '_t' and not self.model.complement)
                or (k[-2:] == '_c' and self.model.complement)):
                setattr(self.model,paramname,params[k])

def LoadEvents(filenames):
    # load all events corresponding to filenames
    events = []
    for fn in filenames:
        try:
            events.append(PoissEvent(fn,'t'))
        except Exception:
            pass
        try:
            events.append(PoissEvent(fn,'c'))
        except Exception:
            pass
    return events