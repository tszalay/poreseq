
class RegionInfo:
    '''Helper class for region data.'''
    
    def __init__(self, region=None):
        '''Parses a region string.
        
        Accepts four formats:
        
        None
        sequence_name
        10000:20000
        sequence_name:10000:20000
        
        and it is up to later processing to fill in unknown fields appropriately.
        '''
        self.start = None
        self.end = None
        self.name = None        

        if region is None:
            return
            
        rs = region.split(':')
        if len(rs) != 2:
            self.name = rs[0]

        if len(rs) > 1:
            self.start = int(rs[-2])
            self.end = int(rs[-1])

class MutationInfo:
    '''Contains information for a single mutation.
    
    Note that orig and mut should be '' and not '.' for insertions/deletions.    
    
    Attributes:
        start(int): 0-based starting base of mutation
        orig(string): sequence of original bases ('ACG')
        mut(string): sequence of mutated bases ('GT')
    '''
    
    def __init__(self, info=None):
        '''Initialize the class, either with blank values, or parsing a string.
        
        If parsing info (eg. reading from file), it takes the format
        <start>     <orig>      <mut>
        using whitespace as a delimiter and ignoring comments. It also coverts
        input '.' strings to ''.
        '''
        self.start = 0
        self.orig = ""
        self.mut = ""
        if info is not None:
            if len(info) == 0 or info[0] == '#':
                self.start = -1
                return
                
            vals = info.split()
            
            if len(vals) != 3:
                self.start = -1
                return
                
            self.start = int(vals[0])
            self.orig = vals[1]
            self.mut = vals[2]
            if self.orig == '.':
                self.orig = ''
            if self.mut == '.':
                self.mut = ''
                
    def __str__(self):
        '''Outputs in a nice format.'''
        original = self.orig
        if len(original) == 0:
            original = '.'
        mutation = self.mut
        if len(mutation) == 0:
            mutation = '.'
            
        return '{}\t{}\t{}'.format(self.start,original,mutation)

class MutationScore:
    '''Contains information for a scored mutation.
    
    Note that orig and mut should be '' and not '.' for insertions/deletions.    
    
    Attributes:
        start(int): 0-based starting base of mutation
        orig(string): sequence of original bases ('ACG')
        mut(string): sequence of mutated bases ('GT')
        score(float): scored change in likelihood if this mutation were made
    '''
    def __init__(self):
        '''Default initializer.'''
        self.start = 0
        self.orig = ""
        self.mut = ""
        self.score = 0
        
    def __str__(self):
        '''Outputs in a nice format.'''
        original = self.orig
        if len(original) == 0:
            original = '.'
        mutation = self.mut
        if len(mutation) == 0:
            mutation = '.'
            
        return '{}\t{}\t{}\t{}'.format(self.start,original,mutation,self.score)
