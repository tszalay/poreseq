import h5py
import pdb
import numpy as np
import copy
from Bio.Seq import Seq



def MakeContiguous(pyobj):
    """Makes all numpy fields of object into memory-contiguous arrays.
    
    Loops through all attributes of pyobj, and makes any numpy.ndarray object
    found contiguous using np.ascontiguousarray (assumes double arrays).
    
    Args:
        pyobj (class): Class to make contiguous
    
    Returns:
        None
        
    Raises:
        None
    """
    attrs = [attr for attr in dir(pyobj) if not callable(attr) and not attr.startswith("__")]
    for attr in attrs:
        val = getattr(pyobj,attr)
        if type(val) == np.ndarray:
            setattr(pyobj,attr,np.ascontiguousarray(val,dtype='f8'))
            
            
def LoadEvents(filenames):
    '''Load all events corresponding to a list of filenames, both t and c'''
    
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
    

class PoissModel:
    """Class containing all trained model parameters for a single event.
    
    Many of the parameters are loaded from the fast5 data, and the remaining
    skip and stay params are provided in the paramfile.
    
    Note that sd_mean and sd_stdv are loaded and scaled, but in the C++
    code they actually get converted to an inverse gaussian distribution.
    
    Attributes:
        level_mean (double x 1024): expected 5mer currents
        level_stdv (double x 1024): width of distribution of level_mean
        sd_mean (double x 1024): expected current noise on 5mers
        sd_stdv (double x 1024): deviation of above
        prob_* (double): probabilities of skips, stays, etc. from params
        name (string): Oxford-given name of model
        complement (boolean): Is this a complement model?
    """
    
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
    """Class containing all data associated with a single read.
    
    In this case, the template and complement of a single-molecule read will
    be stored as separate events.
    
    NOTE: complement reads are flipped by default - that is, the current level
    sequences are reversed, and the model is flipped as well to the reverse
    complement read, so that template and complement reads point in the same
    direction.
        
    Attributes:
        mean (double x N): average of measured current values
        stdv (double x N): deviation of measured current values
        length (double x N): duration of measured current values
        start (double x N): start time of measured current values
        ref_align (double x N): reference base each current level aligns to
        ref_like (double x N): cumulative aligned likelihood at each current level        
        model (PoissModel): model associated with this event
        sequence (string): extracted 2D sequence from fast5
        flipped (boolean): has this event been flipped?
    """
    
    def __init__(self, filename,typ):
        """Loads a template or complement event from a fast5 file.
    
        Loads and scales the event and model from a fast5 file, from one
        side only.
        
        Args:
            filename (string): fast5 filename
            typ (string): 't' or 'c' for template/complement
        
        Returns:
            None
        """  
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
        self.flipped = False
        
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
        '''Calls copy.deepcopy on self and returns it.'''
        
        return copy.deepcopy(self)        
        
    def flip(self, flip_sequence=True):
        """Flips the event and model.
    
        In-place reverses all of the events' double x N fields and maps each model
        state to its reverse complement (ACCGG -> CCGGT).
        
        Args:
            flip_sequence (boolean): whether to reverse-complement self.sequence
        
        Returns:
            None
        """
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
        
        self.flipped = not self.flipped
            
    def mapaligns(self, pairs):
        """Re-maps ref_align using list of pairs of aligned indices between old and new.
    
        Uses the result of swalign or pysam functions to shift the event's
        alignment from one sequence to another, for example:        
            pairs = [(10, 20), (11,21), (12,22), ...]
        specifies that where the current ref_align == 10, it should be shifted
        to 20, 11 -> 21, etc.
            
        Args:
            pairs (list of pairs): ints of the form [(i0,j0),(i1,j1),...]
        
        Returns:
            None
        """
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
        '''Makes numpy arrays in event and model contiguous.'''
        
        MakeContiguous(self)
        MakeContiguous(self.model)
    
    def getrefstats(self):
        """Calculate some simple alignment statistics.
        
        Gets the fraction of current levels that are skips, stays, and insertions.
                
        Args:
            None
        
        Returns:
            Tuple of (skips, stays, insertions)
        """
        
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
        """Uses param dictionary to set events' model params.
    
        Takes a dictionary loaded from a .conf file and uses the values to set
        self.model.prob_*. For example, self.model.prob_skip = params['skip_t']
        (if the event is a template event)
            
        Args:
            params (dict): dictionary of parameters
        
        Returns:
            None
        """
        
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
