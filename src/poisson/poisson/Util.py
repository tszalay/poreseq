
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
