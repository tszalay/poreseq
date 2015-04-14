/*==========================================================
 *
 * AlignData.h
 *
 * Simple struct for holding things required for alignment.
 * Q - But shouldn't this be a class with member functions? 
 *     Since all functions take this as their first argument
 *     anyway.
 * A - Maybe. But because there is no good way to construct
 *     the class given where events/sequence/etc might come
 *     from (Python or Matlab or wherever), this is a better 
 *     representation of the usage.
 *
 *========================================================*/

#ifndef _ALIGNDATA_H_
#define _ALIGNDATA_H_

#include <vector>
#include <map>

#include "Sequence.h"
#include "EventData.h"
#include "AlignUtil.h"

struct AlignData
{
    // the obvious three
    Sequence sequence;
    vector<EventData> events;
    AlignParams params;

    // for holding cached likelihood values that don't change
    map< string, vector<double> > seqlikes;
};

#endif