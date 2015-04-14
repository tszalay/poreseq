#! /usr/bin/env python

import sys
import os
import glob
import h5py

def get_fasta(filename):
    
    # load fast5 file
    f = h5py.File(filename,'r')
    # get the Oxford-called 2D sequence
    seqdata = f['/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq']
    # load sequence (parse from Fastq)
    return seqdata[()].split('\n')[1]
    
def extract_fasta(fast5files, fastafile=None, addpath=False, force=False):

    # did we get a non-empty list?    
    if not fast5files:
        raise Exception('No files specified!')
    
    # if no fasta file specified, find from the files' base path
    if fastafile is None:
        fast5dir = fast5files[0]
        if os.path.isdir(fast5dir):
            fastafile = os.path.normpath(fast5dir)
        else:
            fastafile,_ = os.path.split(fast5dir)
        fastafile += '.fasta'
    
    # check if the file exists
    if os.path.isfile(fastafile) and not force:
        if raw_input('File ' + fastafile + ' exists, overwrite? (y/n): ') != 'y':
            return
        
    fasta = open(fastafile,'w') 
    
    print 'Extracting fasta to ' + fastafile + ' ...     '
    nwrote = 0
    for i,f in enumerate(fast5files):
        try:
            seq = get_fasta(f)
        except:
            continue
        fn = f
        # do we want to put the relative/provided path in the header?
        # this is useful when merging multiple fasta files
        if not addpath:
            _,fn = os.path.split(f)
        fa = '>' + fn + '\n' + seq + '\n'
        fasta.write(fa)
        sys.stdout.write('\r{:3}%'.format(int(100*i/len(fast5files))))
        sys.stdout.flush()
        nwrote = nwrote + 1
        
    print '\rDone, extracted '+ str(nwrote) + ' 2D fasta sequences'
    
    fasta.close()

if __name__ == "__main__":
    
    fast5files = []
    fastafile = None
    addpath = False
    force = False
    
    if len(sys.argv) < 2:
        print "usage: extract_fasta.py [-p] fast5-dir/files [fasta file]"
        print "    (-p option includes fast5 path in fasta header names)"
        print "\nexample: extract_fasta.py data/run_23/*.fast5 out.fasta"
        sys.exit()
        
    for arg in sys.argv[1:]:
        if arg == '-p':
            addpath = True
        elif arg == '-f':
            force = True
        elif arg.find('.fast5') > 0 and os.path.isfile(arg):
            fast5files.append(arg)
        elif arg.find('.fa')>0 or arg.find('.fasta')>0:
            fastafile = arg
        elif os.path.isdir(arg):
            # glob the fast5 files
            fast5dir = os.path.join(arg,'*.fast5')
            fast5files += glob.glob(fast5dir)
        
    extract_fasta(fast5files,fastafile,addpath,force)
