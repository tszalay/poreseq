/*==========================================================
 *
 * MakeMutations.cpp
 *
 * Functions for taking a bunch of strands of data, a sequence,
 * and an optional list of mutations, and optimizes the sequence
 * by testing mutations iteratively.
 *
 *========================================================*/


#include "Mutations.h"
#include "Alignment.h"

// define comparison operator here, as one might think
bool operator< (const MutScore& ms1, const MutScore& ms2)
{ return ms1.score > ms2.score; }

// calculate the scores for the given list of mutations, using all
// of the events we have. does it event-by-event instead of mut-by-mut,
// so that if one wanted to parallelize it by doing the different events
// simultaneously, it would not be difficult (also saves memory)
vector<MutScore> ScoreMutations(AlignData& data, const vector<MutInfo>& muts)
{
    // initialize scored output vector as copy of MutInfos
    vector<MutScore> mutscores = vector<MutScore>(muts.begin(),muts.end());
    
    if (data.params.verbose)
    {
        cerr << "Scoring (" << data.params.scoring_width << ")";
        cerr.flush();
    }

    vector<Alignment> alignments;
    for (int i=0; i<data.events.size(); i++)
        alignments.push_back(Alignment(data.sequence, data.events[i], data.params));
    
    for (int j=0; j<alignments.size(); j++)
    {
        // now fill/refill only this alignment, and add the scores from it
        alignments[j].update(data.sequence);
        
        for (int i=0; i<mutscores.size(); i++)
        {
            // sanity check to make sure mutation isn't crappy
            if (muts[i].start > data.sequence.bases.size())
                continue;
            // create mutated sequence
            Sequence mutseq(data.sequence,muts[i]);
            // and start the mutated alignments
            mutscores[i].score += alignments[j].scoreMutation(muts[i],mutseq);
        }
        // and now clear this alignment, to save memory
        alignments[j].clear();
        if (data.params.verbose)
        {
            cerr << ".";
            cerr.flush();
        }
    }

	if (data.params.verbose)
	{
		cerr << endl;
		cerr.flush();
    }
	
    return mutscores;
}

// applies the given mutations in order of scores, but only applying ones
// that are far enough apart that they're unlikely to affect one another;
// ones that are too close together, it re-scores and re-applies
int MakeMutations(AlignData& data, vector<MutScore> muts)
{
    // how far from a made mutation do we assume scores changed etc?
    const int mutspc = 10;
    
	// how many bases did we mutate?
	int mutbases = 0;
    
    // and sort them by score, descending
    sort(muts.begin(), muts.end());
    // remove all bad (negative) mutations
    while (muts.size() > 0 && muts.back().score < 0)
        muts.pop_back();
    // done if no mutations left
    if (muts.size() == 0)
        return 0;

    if (data.params.verbose)
    {
        cerr << "Testing " << muts.size() << " mutations..." << endl;
        cerr.flush();
    }
    
    // vector in which to put ones we didn't use
    vector<MutInfo> mutextra;
    
    // and run through them in order, making them in the sequence
    for (int i=0; i<muts.size(); i++)
    {
        // was it previously invalidated?
        if (muts[i].score < 0)
        {
            mutextra.push_back(MutInfo(muts[i]));
            continue;
        }
        // otherwise, make the mutation
        data.sequence = Sequence(data.sequence,muts[i]);
        
        if (data.params.verbose > 1)
        {
            cerr << "Kept mutation " << i << " at " << muts[i].start 
                    << " of " << muts[i].orig.size() << " to " << muts[i].mut.size() 
                    << " with score " << muts[i].score << endl;
            cerr.flush();
        }
        // add how many bases we mutated
		mutbases += max(muts[i].orig.size(),muts[i].mut.size());

        // we kept a mutation, so we need to modify all subsequent 
        // mutations to account for change
        for (int j=i+1; j<muts.size(); j++)
        {
            // does it overlap with padded, then run it through on the next round instead
            int minind = max(muts[i].start,muts[j].start);
            int maxind = min(muts[i].start+muts[i].mut.size(),muts[j].start+muts[j].mut.size());
            if (minind < maxind+mutspc && muts[j].score > 0)
            {
                // overlaps, do it later
                muts[j].score = -1;
                continue;
            }

            // doesn't overlap, if it starts after, offset the start
            if (muts[j].start >= muts[i].start+muts[i].orig.size())
                muts[j].start += (muts[i].mut.size() - muts[i].orig.size());
        }
    }
    
    if (mutextra.size() > 10)
        mutbases += MakeMutations(data,ScoreMutations(data,mutextra));
	
	return mutbases;
}

vector<double> ScoreAlignments(AlignData& data, double* likes)
{
    // loop through and do each alignment
    // this in-place updates each event and then
    // clears the alignment to save memory
    vector<double> scores;

    vector<Alignment> alignments;
    for (int i=0; i<data.events.size(); i++)
        alignments.push_back(Alignment(data.sequence, data.events[i], data.params));
    
    for (int i=0; i<alignments.size(); i++)
    {
        // calculate the matrix
        alignments[i].fillColumns();
        // and save the backtrace in the event, getting its pointer
        EventData* event = alignments[i].backtrace();
        // as well as the output score
        scores.push_back(alignments[i].getMax());
        // and the likelihood, if requested
        if (likes > 0)
        {
            // save the aligned likelihoods into lpr
            double lastlik = 0;
            int refind = 1;
            // loop through and find ref_aligns
            for (int j=0; j<event->length; j++)
            {
                // save these refinds
                if (event->ref_align[j] > 0)
                {
                    // if ref_align[j] == refind, it just updates lastlik
                    for (int k=refind; k<event->ref_align[j]; k++)
                        likes[k+1] += lastlik;
                    lastlik = event->ref_like[j];
                    refind = event->ref_align[j];
                }
            }
            // and add the last ones on
            for (int k=refind; k<data.sequence.states.size()+3; k++)
                likes[k+1] += lastlik;
        }
        
        alignments[i].clear();
    }
    
    return scores;
}