/*==========================================================
 *
 * Mutations.h
 *
 * Header file for all mutation and alignment-related functions.
 *
 *========================================================*/

#ifndef _MUTATIONS_H_
#define _MUTATIONS_H_

#include <vector>

#include "Sequence.h"
#include "AlignData.h"

// in FindMutations.cpp
vector<MutInfo>     FindMutations(AlignData& data, const vector<Sequence>& seqs);
vector<MutInfo>     FindPointMutations(AlignData& data);

// in MakeMutations.cpp
int                 MakeMutations(AlignData& data, vector<MutScore> muts);
vector<MutScore>    ScoreMutations(AlignData& data, const vector<MutInfo>& muts);
vector<double>      ScoreAlignments(AlignData& data, double* likes);

#endif