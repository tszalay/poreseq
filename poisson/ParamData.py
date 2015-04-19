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
            f.write('{} = {}\n'.format(p,params[p]))
            
def VaryParams(params):
    '''Return a list of param dicts with modifications'''

    paramlist = []
    fac = 1.2
    
    for p in params:
        # only vary relevant (non-configurational) params
        if p[-2:]=='_t' or p[-2:]=='_c':
            # multiply and divide it by a certain factor
            newparams = params.copy()
            newparams[p] = params[p]*fac
            paramlist.append(newparams)
            newparams = params.copy()
            newparams[p] = params[p]/fac
            paramlist.append(newparams)
            
    return paramlist