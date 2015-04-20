#! /usr/bin/env python

import sys
import os
import h5py

def get_fasta(filename):
    
    # load fast5 file
    f = h5py.File(filename,'r')
    # get the Oxford-called 2D sequence
    seqdata = f['/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq']
    seq = seqdata[()].split('\n')[1]
    # and close the h5py file, just in case
    f.close()
    return seq
    
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
        sys.stderr.write('File exists, skipping...\n')
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
