/*==========================================================
 *
 * EventUtil.h
 *
 * Useful functions for event alignment.
 *
 *========================================================*/

#ifndef _EVENTUTIL_H_
#define _EVENTUTIL_H_

#include "AlignData.h"
#include "EventData.h"
#include "Sequence.h"
#include "swlib.h"

SWAlignment MapAlignments(AlignData& data, const Sequence& newseq);

//void MapAlignment(EventData& event, const SWAlignment& align);
/*void mapaligns(vector<EventData>& events, const SWAlignment& align);
void realign(Sequence& seq, EventData& event, AlignParams params);
void seedaligns(Sequence& seq, vector<EventData>& events, AlignParams params);
*/

#endif