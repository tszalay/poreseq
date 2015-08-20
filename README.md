# PoreSeq

The PoreSeq utility is an open source program and Python library for de novo sequencing, consensus and variant calling on data from Oxford Nanopore Technologies' MinION platform. Features include:

* de novo error correction without reference using overlap alignment
* reference error correction
* scoring known sequence variants on a given dataset
* straightforward subdivision of processing for cluster/parallel tasks

## Install
The recommended way to install poreseq is from the source using git and pip:

1. Checkout the source: `git clone git://github.com/tszalay/poreseq.git`
2. In the main directory, run `pip install -e . --user` or `sudo pip install -e .`

Note that for a userspace install, `~/.local/bin` needs to be on your PATH.
The Python requirements are handled by setuptools, while some of the utility scripts depend on the LAST aligner and the Celera/PBcR assembler. If you would like
to change these requirements, edit the scripts in the scripts/ folder (and reinstall). Note that the poreseq_assemble script, in addition to
needing PBcR, also needs the correct path to poreseq.spec and should be updated as such.
   
## Examples
The basic aspects of the usage are described below. For any command, the -h flag displays help.

Sequence extraction from fast5 files:

* `poreseq extract /media/run-33 ./allreads.fasta`

Overlap alignment via LAST (outputs to allaligns.bam):

* `poreseq_align ./allreads.fasta ./allreads.fasta allaligns`

Split fasta into 10kb regions in 16 files for processing:

* `poreseq split ./allreads.fasta -n 16 -R 10000`

Error correction (correct the regions in the 15th file as defined by the splitting above):

* `poreseq consensus -p params.conf ./allreads.fasta ./allaligns.bam /media/run-33 -R allreads.44.region -o corr.44.fasta`

Error correction (manually correct a specific region of a read from allreads):

* `poreseq consensus -p params.conf ./allreads.fasta ./allaligns.bam /media/run-33 -r ch_23_read_11.fast5:1000:2000 -o corr.fasta`

Merge previously split and corrected sequences:

* `poreseq merge merged.fasta corr.*.fasta`

Assemble sequences using Celera assembler (with 4 threads):

* `poreseq_assemble merged.fasta genomeSize=48500 -t 4`

Train skip and stay parameters on bases 46kb-48kb of a known reference (with 16 threads):

* `poreseq train reference.fasta alignment.bam reads.fasta /media/run-data -n 16 -i 30 -p starting_params.conf -r 46000:48000`

Score mutations listed in mutations.txt against reference:

* `poreseq variant reference.fasta alignment.bam /media/run-data -m muts.txt -p params.conf`

Example mutations.txt (start, original, mutated):

```
11      AC    G
1994    .     C
2301    A     .
```

(note that start index is 0-based)

Example of Python package usage:
```
import poreseq

params = poreseq.LoadParams('./params.conf')
reginfo = poreseq.RegionInfo('10000:20000')  # poreseq.RegionInfo() for whole strand
pa = poreseq.LoadAlignedEvents('reference.fasta','alignment.bam','/media/run-data',reginfo,params)
# pa is now a PSAlign class, pa.sequence is the consensus/ref sequence
# pa.events are the loaded events and pa.params is the same params as loaded above
# get likelihood score of each event, also realign events to ref.
pa.ScoreEvents()
# display which reference bases are aligned to by current levels 2000:3000 for event 2
# note that these are 1-based for internal reasons
print pa.events[2].ref_align[2000:3000]
```

At the time of writing, “poreseq consensus” takes around 2 minutes to error-correct a 1 kb region at 10X coverage, meaning that it can take tens of hours on a single CPU to error-correct the fragments before assembling λ DNA. If the coverage is sufficient, it is more efficient to assemble without doing poreseq error correction, and then using poreseq to align and refine the resulting assembly as described above. Additionally, if only a specific region of the genome is desired, poreseq can simply error correct that region, greatly reducing the time required.

