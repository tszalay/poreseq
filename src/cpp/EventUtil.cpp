/*==========================================================
 *
 * EventUtil.cpp
 *
 * A few handy functions for events and alignments.
 *
 *========================================================*/

#include "EventUtil.h"
#include "Alignment.h"

SWAlignment MapAlignments(AlignData& data, const Sequence& newseq)
{
    // align the two sequences, fill in insertions/deletions with the
    // previous index
    SWAlignment align = fillinds(swfull(data.sequence.bases,newseq.bases));
    data.sequence = newseq;

    const vector<int>& inds1 = align.inds1;
    const vector<int>& inds2 = align.inds2;
    
    for (int i=0; i<data.events.size(); i++)
    {
        EventData& event = data.events[i];
        // now, we want to re-populate the ref_align with new ref_align
        // values pointing to inds2, this will be approximate but sufficient
        // -----------------
        // in other words, we want to evaluate inds2 = f(inds1) at the points
        // from original refal
        for (int j=0; j<event.ref_align.size(); j++)
        {
            int refal = (int)event.ref_align[j];
            if (refal < inds1.front() || refal > inds1.back())
            {
                // no data, set to 0
                event.ref_align[j] = 0;
                continue;
            }
            // otherwise, find nearest value in inds1 (first value not less than)
            auto it = std::lower_bound(inds1.begin(),inds1.end(),refal);
            // get the index
            int ind = it - inds1.begin();
            // and set that inds2 value, if found
            if (ind < inds2.size())
                event.ref_align[j] = inds2[ind];
            else
                event.ref_align[j] = 0;
        }
    
        // and update the internal data, pretty key
        event.updaterefs();
    }

    return align;
}
