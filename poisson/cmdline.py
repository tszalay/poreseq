from Mutate import Mutate
from Util import *
from ParamData import *
from extract_fasta import extract_fasta
from split_fasta import split_fasta, split_regions
from merge_fasta import merge_fasta

import string
import argparse
import sys
import os
import stat
import glob

import numpy as np
from multiprocessing import Pool
from Bio import SeqIO

def main():

    parser = argparse.ArgumentParser(prog='poisson')
    subparsers = parser.add_subparsers(help='Nanopore sequence consensus tool')
    
    # create consensus-calling parser
    parse_cons = subparsers.add_parser('consensus', help='run consensus algorithm using alignment')
    parse_cons.add_argument('ref', help='reference fasta file')
    parse_cons.add_argument('bam', help='input BAM file')
    parse_cons.add_argument('dir', help='root fast5 directory')
    parse_cons.add_argument('-r', '--regions', default=None, nargs='+',
                        help='regions to correct (eg. 1000:3000 or header_name:1000:3000)')
    parse_cons.add_argument('-R', '--region-file', default=None, nargs='+',
                        help='file(s) containing region strings, one per line')
    parse_cons.add_argument('-p', '--params', default=None,
                        help='parameter file to use')
    parse_cons.add_argument('-v', '--verbose', action="count", default=0,
                        help='output verbosity (0-2)')
    parse_cons.add_argument('-o', '--output', default=None,
                        help='output fasta file ("-" for stdout)')
    parse_cons.add_argument('-T', '--test', action="store_true", default=False,
                        help='test mode: seed with loaded sequence, output score as well')
    parse_cons.set_defaults(func=consensus)
    
    # and the variant-calling parser
    parse_var = subparsers.add_parser('variant', help='call sequence variants')
    parse_var.add_argument('ref', help='reference fasta file')
    parse_var.add_argument('bam', help='input BAM file')
    parse_var.add_argument('dir', help='root fast5 directory')
    parse_var.add_argument('vars', help='fasta of variant sequences to test')
    
    parse_var.add_argument('-r', '--regions', default=None, nargs='+',
                        help='regions to correct (eg. 1000:3000 or header_name:1000:3000)')
    parse_var.add_argument('-p', '--params', default=None,
                        help='parameter file to use')
    parse_var.add_argument('-v', '--verbose', action="count", default=0,
                        help='output verbosity (0-2)')
    parse_var.set_defaults(func=variant)
                        
    # and the training parser
    parse_train = subparsers.add_parser('train', help='train model parameters on data')
    parse_train.add_argument('ref', help='reference fasta file')
    parse_train.add_argument('bam', help='input BAM file')
    parse_train.add_argument('dir', help='root fast5 directory')

    parse_train.add_argument('-i', '--iter', type=int, default=30,
                        help='number of training iterations')
    parse_train.add_argument('-n', '--threads', type=int, default=4,
                        help='number of threads to use')
    parse_train.add_argument('-p', '--params', default=None,
                        help='initial parameter file to use')
    parse_train.add_argument('-r', '--region', default=None, 
                        help='region to train on (eg. 1000:3000 or header_name:1000:3000)')
    parse_train.set_defaults(func=train)
    
    # short utility parsers: split
    parse_split = subparsers.add_parser('split', help='split fasta files into chunks')
    parse_split.add_argument('fasta', help='fasta file')
    parse_split.add_argument('-R', '--region-length', type=int, default=None,
                        help='output region strings split into this many bases')
    group = parse_split.add_mutually_exclusive_group(required=True)
    group.add_argument('-n', '--num-files', type=int, default=None,
                        help='number of files to split into')
    group.add_argument('-m', '--per-file', type=int, default=None,
                        help='number of sequences per file')
    parse_split.set_defaults(func=split)

    # short utility parsers: merge
    parse_merge = subparsers.add_parser('merge', help='merge corrected fasta files')
    parse_merge.add_argument('fasta_out', help='output fasta filename')
    parse_merge.add_argument('fasta_in', nargs='+', help='fasta files to merge')
    parse_merge.set_defaults(func=merge)

    # and extract
    parse_ext = subparsers.add_parser('extract', help='extract fasta from fast5')
    parse_ext.add_argument('dirs', help='fast5 directories', nargs='+')
    parse_ext.add_argument('fasta', help='output fasta')
    parse_ext.add_argument('-p', '--path', action="store_true", default=False,
                            help='use rel. path as fasta header (instead of just filename)')

    parse_ext.set_defaults(func=extract)


    
    args = parser.parse_args()
    args.func(args)
        
        
def consensus(args):

    # load the parameters file first
    params = LoadParams(args.params)


    # if we specified region files, load regions from there, directly, no modifications
    if args.region_file is not None:
        if args.regions is None:
            args.regions = []
        for rf in args.region_file:
            if os.path.isfile(rf):
                args.regions += [x.strip() for x in open(rf).readlines()]

    # open reference sequence and see how many there are
    # (if we gave a file with multiple sequences and no region string, generate one)
    # and we should make sure to split the regions up, if they weren't directly given
    refs = SeqIO.index(args.ref, "fasta")
    if args.regions is None:
        args.regions = []
        for refid in refs:
            # load the reference sequence
            refseq = refs[refid]
            # check if we need to take sub-slices
            if 'max_length' in params and len(refseq) > params['max_length']:
                dl = int(params['max_length']) - 1000
                istart = 0
                iend =int( params['max_length'])
                # step through max_length-1000 at a time (leave nice overlap)
                while istart < iend:
                    args.regions.append('{}:{}:{}'.format(refid,istart,iend))
                    iend = min(iend+dl,len(refseq))
                    istart = min(istart+dl,len(refseq))
            else:
                # or just append the name directly
                args.regions.append(refid)
                
    # now create output file for writing (or stdout)
    if args.output is None:
        fastabase = os.path.splitext(args.ref)[0]
        outfile = open(fastabase+'.corr.fasta','w')
    elif args.output == '-':
        outfile = sys.stdout
    else:
        outfile = open(args.output,'w')
            
    # now loop through and refine sequences
    # (continue on errors)
    for region in args.regions:
        try:
            seq = Mutate(args.ref,args.bam,args.dir,params=params,region=region,
                         test=args.test,verbose=args.verbose)
        except Exception as e:
             sys.stderr.write('Skipping {}: {}\n'.format(region,str(e)))
             continue
                         
        # test mode output returns accuracy as well
        if args.test:
            (seq, acc) = seq
            
        # print corrected fasta string with region as header
        if args.test:
            region += ' [' + str(round(acc,2)) + ']'
            
        outfile.write('>{}\n{}\n'.format(region,seq))
        outfile.flush()
        

        
def variant(args):
    pdb.set_trace()
    pass
    

# nice trick from stackoverflow to allow function pickling
class trainhelper(object):
    def __init__(self, _args):
        self.args = _args
    def __call__(self, params):
        return Mutate(self.args.ref,self.args.bam,self.args.dir,params=params,region=self.args.region,
               test=True,verbose=False)

def train(args):
    
    # load the specified parameter file
    params = LoadParams(args.params)
    
    # and start the iterations
    for i in range(args.iter):
        # first, generate modified parameters
        paramlist = VaryParams(params)
        
        # use a python multiprocess pool to run mutation on all of them
        pool = Pool(processes=args.threads)
        seqs = pool.map(trainhelper(args), paramlist)
        accs = [s[1] for s in seqs]
        
        # and update the global params
        params = paramlist[np.argmax(accs)]
        # and intermittently save it
        SaveParams('train_best.conf',params)
        # and output to stderr with an update
        sys.stderr.write('Best at iter {}: {}\n'.format(i+1,max(accs)))

def extract(args):
    
    fast5files = []
    for d in args.dirs:
        fast5files += glob.glob(os.path.join(d,'*.fast5'))
    extract_fasta(fast5files,args.fasta,args.path,False)

    
def split(args):
    if args.region_length is None:
        split_fasta(args.fasta,args.num_files,args.per_file)
    else:
        split_regions(args.fasta,args.region_length,args.num_files,args.per_file)

def merge(args):

    merge_fasta(args.fasta_in,args.fasta_out)
