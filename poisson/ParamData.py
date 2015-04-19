def LoadParams(filename):
    '''Load parameter configuration file'''
    params = {}
    if filename is None:
        return params
        
    with open(filename) as f:
        splitlines = [l.split('=') for l in f.readlines()]
        for sl in splitlines:
            if len(sl) == 2:
                pname = sl[0].strip()
                params[pname] = float(sl[1])
                
    return params
    
def SaveParams(filename, params):
    '''Save parameter configuration file'''
    with open(filename,'w') as f:
        for p in params:
            f.write('{} = {}\n',p,params[p])