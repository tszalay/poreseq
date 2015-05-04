
class RegionInfo:

    def __init__(self, region):

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
    def __init__(self, info=None):
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
        original = self.orig
        if len(original) == 0:
            original = '.'
        mutation = self.mut
        if len(mutation) == 0:
            mutation = '.'
            
        return '{}\t{}\t{}'.format(self.start,original,mutation)

class MutationScore:
    def __init__(self):
        self.start = 0
        self.orig = ""
        self.mut = ""
        self.score = 0
        
    def __str__(self):
        original = self.orig
        if len(original) == 0:
            original = '.'
        mutation = self.mut
        if len(mutation) == 0:
            mutation = '.'
            
        return '{}\t{}\t{}\t{}'.format(self.start,original,mutation,self.score)
