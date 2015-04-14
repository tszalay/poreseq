/*==========================================================
 *
 * FindMutations.cpp
 *
 * Find likely mutations given sequences, seq, events, params
 * where the mutations are drawn from the sequences array
 *
 *========================================================*/

#include "Mutations.h"
#include "EventUtil.h"

#include <algorithm>

using namespace std;

// simple arg_max function for vectors
int max_index(vector<double> vec)
{
    return std::max_element(vec.begin(),vec.end())-vec.begin();
}

// finds likely mutations using seqs as seed strands for said mutations
vector<MutInfo> FindMutations(AlignData& data, const vector<Sequence>& seqs)
{
    // start by re-aligning events as default alignment
    // (and keeping their aligned likelihoods)
    vector<double> seqreflike(data.sequence.bases.size(),0);
    ScoreAlignments(data, seqreflike.data());
    
    vector< vector<double> > alllikes;
    vector<SWAlignment> seqals;

    if (data.params.verbose)
        cerr << "Finding mutations";
    
    for (int i=0; i<seqs.size(); i++)
    {
        // now remap the alignments for the new events
        AlignData newdata(data);
        SWAlignment align = MapAlignments(newdata, seqs[i]);
        // and realign using those alignments, getting new likes in the process
        // use the cache in data to store the values
        vector<double>& reflikes = data.seqlikes[seqs[i].bases];
        if (reflikes.size() == 0)
        {
            reflikes = vector<double>(seqs[i].bases.size(),0);
            ScoreAlignments(newdata,reflikes.data());
        }
        
        // decrement inds1 and inds2, because i did it in matlab code
        for (int j=0; j<align.inds1.size(); j++)
        {
            align.inds1[j]-=2;
            align.inds2[j]-=2;
        }
        
        // and remove invalid indices
        while (align.inds1[0] < 0 || align.inds2[0] < 0)
        {
            align.inds1.erase(align.inds1.begin());
            align.inds2.erase(align.inds2.begin());
        }
        // now everything in align.inds1,inds2 should be valid rel.
        // to the sequences (in the range 0, n-1 for that sequence)
        vector<double> alref1;
        vector<double> alref2;
        // set these guys to the swaligned-likelihood values
        for (int j=0; j<align.inds1.size(); j++)
        {
            alref1.push_back(seqreflike[align.inds1[j]]);
            alref2.push_back(reflikes[align.inds2[j]]);
        }
        // take their difference
        for (int j=alref1.size()-1; j>0; j--)
        {
            alref1[j] -= alref1[j-1];
            alref2[j] -= alref2[j-1];
        }
        alref1[0] = 0;
        alref2[0] = 0;
        // now do CUSUM, re-integrate, but clamping negatives to 0
        vector<double> dlikes(alref1.size());
        double cusum = 0;
        for (int j=0; j<alref1.size(); j++)
        {
            double d = alref2[j] - alref1[j];
            cusum += d;
            if (cusum < 0) cusum = 0;
            dlikes[j] = cusum;
            // also set to 0 if they are basically the same
            // this means the alignment there is the same
            if (abs(alref1[j]-alref2[j]) < 1e-5) dlikes[j] = 0;
        }
        
        // save all the stuff for mutation finding
        alllikes.push_back(dlikes);
        seqals.push_back(align);

        if (data.params.verbose)
        {
            cerr << ".";
            cerr.flush();
        }

    }

    if (data.params.verbose)
        cerr << endl;
        
    // now continue on and actually find the mutations
    vector<MutInfo> mutations;

    // limit the number we return to aroudn here, arbitrarily
    while (mutations.size() < data.sequence.bases.size()/3)
    {
        // find maximum likelihood out of all dlikes
        vector<double> lmax(alllikes.size(),0);
        for (int i=0; i<alllikes.size(); i++)
            lmax[i] = alllikes[i][max_index(alllikes[i])];

        // find the indices where it occured
        int imax = max_index(lmax);
        int ind = max_index(alllikes[imax]);
        
        vector<double>& dlike = alllikes[imax];
        
        // if it's a pretty small value, don't bother saving any more
        if (dlike[ind] < 0.25)
            break;
        
        // now find previous and next zero
        // iterator to max elt.
        auto maxiter = dlike.begin()+ind;
        // next zero after max element @ i1
        int i1 = std::find(maxiter,dlike.end(),0)-dlike.begin();
        // go in reverse to find previous 0 @ i0
        int i0 = std::find(reverse_iterator<vector<double>::iterator>(maxiter+1),dlike.rend(),0).base()-dlike.begin();
        i0--; // subtract one b/c
        
        // set to beginning or end if not found
        if (i0 < 0) i0 = 0;
        if (i1 < 0) i1 = 0;
        if (i0 >= dlike.size()) i0 = dlike.size()-1;
        if (i1 >= dlike.size()) i1 = dlike.size()-1;
        
        // now pull out corresponding indices in the two sequences
        // these are valid based on their original values
        int start1 = seqals[imax].inds1[i0];
        int start2 = seqals[imax].inds2[i0];
        int end1 = seqals[imax].inds1[ind];
        int end2 = seqals[imax].inds2[ind];
        
        // now actually extract the mutation info
        MutInfo mut;
        mut.start = start1;
        // save original and mutated strings
        mut.orig = data.sequence.bases.substr(start1,end1-start1);
        mut.mut = seqs[imax].bases.substr(start2,end2-start2);
        // and trim unmutated bases from front (increase start)
        while (mut.orig.size() > 0 && mut.mut.size() > 0
                && *mut.orig.begin() == *mut.mut.begin())
        {
            mut.orig.erase(mut.orig.begin());
            mut.mut.erase(mut.mut.begin());
            mut.start++;
        }
        // and from the back (start stays the same)
        while (mut.orig.size() > 0 && mut.mut.size() > 0
                && *mut.orig.rbegin() == *mut.mut.rbegin())
        {
            mut.orig.erase(mut.orig.end()-1);
            mut.mut.erase(mut.mut.end()-1);
        }
        
        // only do something if we're mutating something
        // if they were the same string, they're empty from above
        if (mut.orig.size() > 0 || mut.mut.size() > 0)
            mutations.push_back(mut);
        
        // now clear that range of dlike
        std::fill(dlike.begin() + i0, dlike.begin() + i1 + 1, 0);
    }
    
    return mutations;
}


// returns a list of the scores of all possible point mutations
// (including probably lots of duplicates)
vector<MutInfo> FindPointMutations(AlignData& data)
{
    // run through and test all insertions/mutations/deletions and score them

    // symbols for insertions/mutations, make static const for persistence
    static char* bases = "ACGT";
    
    vector<MutInfo> muts;
    
    for (int i=0; i<data.sequence.states.size(); i++)
    {
        MutInfo mut;
        mut.start = i;
        
        // deletion
        mut.orig = data.sequence.bases[i];
        mut.mut = "";
        muts.push_back(mut);
        
        // mutations
        mut.orig = data.sequence.bases[i];
        for (int j=0; j<4; j++)
        {
            // skip over the non-mutation
            if (data.sequence.bases[i] == bases[j])
                continue;
            mut.mut = bases[j];
            muts.push_back(mut);
        }
        
        // insertions
        mut.orig = "";
        for (int j=0; j<4; j++)
        {
            mut.mut = bases[j];
            muts.push_back(mut);
        }
    }
    
    if (data.params.verbose)
        cerr << "Point ";
    
    return muts;
}
