def LoadParams(filename):
    '''Load paramer configuration file'''
    params = {}
    with open(filename) as f:
        splitlines = [l.split('=') for l in f.readlines()]
        for sl in splitlines:
            if len(sl) == 2:
                pname = sl[0].strip()
                params[pname] = float(sl[1])
                
    return params