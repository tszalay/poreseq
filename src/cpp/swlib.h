/*==========================================================
 *
 * swlib.h
 *
 * Functions for computing a striped Smith-Waterman alignment
 * when the rough alignment is known through other means
 * (or just for computing using full matrix otherwise)
 *
 *========================================================*/

#ifndef _SWLIB_H_
#define _SWLIB_H_

#include <iostream>
#include <vector>
#include <stdint.h>
#include <cmath>

using namespace std;

const int score_match = 5;
const int score_mismatch = -4;
const int score_insert = -8;

// struct for holding results
struct SWAlignment
{
    int score;
    double accuracy;
    
    vector<int> inds1;
    vector<int> inds2;
};

// full matrix alignment
SWAlignment swfull(const string& seq1, const string& seq2);
// striped alignment
SWAlignment swfast(const string& seq1, const string& seq2, double al_m, double al_b, int width);
// fill in 0-indices with correct indices
SWAlignment fillinds(SWAlignment align);

#endif