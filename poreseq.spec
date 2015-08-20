# based on Loman's oxford.spec with the new PBcR

#doFragmentCorrection = 0
#assembleCoverage = 5
#assembleMinCoverage = 2

# decrease mer size
merSize=14

# use falcon_sense consensus with adjusted parameters for lower quality reads
falconForce=1
falconOptions=--max_n_read 200 --min_idt 0.50 --output_multi --local_match_count_threshold 0

# adjust assembly parameters to overlap at high error rate since the corrected reads are not 99% like pacbio
asmOvlErrorRate = 0.3
asmUtgErrorRate = 0.3
asmCgwErrorRate = 0.3
asmCnsErrorRate = 0.3
asmOBT=0
batOptions=-RS -CS
utgGraphErrorRate = 0.3
utgMergeErrorRate = 0.3
