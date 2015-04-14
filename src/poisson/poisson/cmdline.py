import argparse
from DoMutation import *
from Util import *

def main():
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
                        help='Test mode, seed with loaded sequence')

    args = parser.parse_args()
    
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
        # trim ends as requested
        if len(seq) > 2*args.t:
            seq = seq[args.t:-args.t]
            
        # print corrected fasta string with region as header
        sys.stdout.write('>{}\n{}\n'.format(region,seq))

    