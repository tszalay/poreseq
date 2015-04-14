import argparse
from DoMutation import *
from Util import *
from ParamData import *
import sys

def main():

    parser = argparse.ArgumentParser(prog='poisson')
    parser.add_argument('mode', help='command to execute')
    # set options here!
    
    args = parser.parse_args()
    
    if args.mode == 'consensus':
        consensus(sys.argv[2:end])
    elif args.mode == 'variant':
        variant(sys.argv[2:end])
    elif args.mode == 'train':
        train(sys.argv[2:end])
    else:
        sys.stderr.write('Invalid mode: {}\n'.format(args.mode))

    sys.exit()
    

        
 def consensus(argv):
 
    parser = argparse.ArgumentParser(prog='poisson')
    parser.add_argument('ref', help='reference fasta file')
    parser.add_argument('bam', help='input BAM file')
    parser.add_argument('dir', help='root fast5 directory')
    parser.add_argument('-r', '--regions', default=None, nargs='+',
                        help='regions to correct (eg. 1000:3000 or header_name:1000:3000)')
    parser.add_argument('-m', '--overlap', type=int, default=300,
                        help='minimum overlap of aligned reads, in bases')
    parser.add_argument('-c', '--coverage', type=int, default=25,
                        help='maximum coverage depth for aligned reads')
    parser.add_argument('-t', type=int, default=25,
                        help='number of bases to trim from ends')
    parser.add_argument('-p', '--params', default=None,
                        help='parameter file to use')
    parser.add_argument('-v', '--verbose', action="count", default=0,
                        help='output verbosity (0-2)')
    parser.add_argument('-T', '--test', action="store_true", default=False,
                        help='test mode: seed with loaded sequence, output score as well')
                        
    args = parser.parse_args(argv)
    
    # open reference sequence and see how many there are
    # (if we gave a file with multiple sequences and no region string, generate one)
    refs = SeqIO.index(args.ref, "fasta")
    if args.regions is None:
        args.regions = list(refs.keys())

    # now loop through and refine sequences
    for region in args.regions:
        seq = DoMutation(args.ref,args.bam,args.dir,paramfile=args.params,region=region,
                         overlap=args.overlap,maxcoverage=args.coverage,test=args.test,
                         verbose=args.verbose)
                         
        # test mode output returns accuracy as well
        if test:
            (seq, acc) = seq
            
        # trim ends as requested
        if len(seq) > 2*args.t:
            seq = seq[args.t:-args.t]
            
        # print corrected fasta string with region as header
        if args.test:
            region += ' [' + str(round(acc,2)) + ']'
            
        sys.stdout.write('>{}\n{}\n'.format(region,seq))

        
def variant(argv):

    parser = argparse.ArgumentParser(prog='poisson')
    parser.add_argument('ref', help='reference fasta file')
    parser.add_argument('bam', help='input BAM file')
    parser.add_argument('dir', help='root fast5 directory')
    parser.add_argument('vars', help='variant sequences to test')
    
    parser.add_argument('-r', '--regions', default=None, nargs='+',
                        help='regions to correct (eg. 1000:3000 or header_name:1000:3000)')
    parser.add_argument('-m', '--overlap', type=int, default=300,
                        help='minimum overlap of aligned reads, in bases')
    parser.add_argument('-c', '--coverage', type=int, default=25,
                        help='maximum coverage depth for aligned reads')
    parser.add_argument('-p', '--params', default=None,
                        help='parameter file to use')
    parser.add_argument('-v', '--verbose', action="count", default=0,
                        help='output verbosity (0-2)')
                        
    args = parser.parse_args(argv)

    
    
def train(argv):
    
    parser = argparse.ArgumentParser(prog='poisson')
    
    # use the arguments exactly as given, and then call them on each parameter file
    parser.add_argument('ref', help='reference fasta file')
    parser.add_argument('bam', help='input BAM file')
    parser.add_argument('dir', help='root fast5 directory')

    parser.add_argument('-i', '--iter', type=int, default=30,
                        help='number of training iterations')
    parser.add_argument('-n', '--threads', type=int, default=4,
                        help='number of threads to use')
    parser.add_argument('-m', '--overlap', type=int, default=300,
                        help='minimum overlap of aligned reads, in bases')
    parser.add_argument('-c', '--coverage', type=int, default=25,
                        help='maximum coverage depth for aligned reads')
    parser.add_argument('-p', '--params', default=None,
                        help='parameter file to use')

    args = parser.parse_args(argv)
    
    # now load the specified parameter file
    params = LoadParams(args.params)
    
    # and start the iterations
    for i in range(args.iter):
        # first, generate modified parameter files and save them
        
        # next, use parallel to spin up all of them at once
        
        # finally, pull out the accuracy from the results
        
        # and load corresponding params file
        
        
        
        