/*==========================================================
 *
 * swlib.cpp
 *
 * Algorithm for computing a striped Smith-Waterman alignment
 * when the rough alignment is known through other means
 * (or just for computing using full matrix otherwise)
 *
 * Arguments: seq0, seq1, paired inds, stripe width
 *
 *========================================================*/

#include "swlib.h"

#include <algorithm>

using namespace std;

SWAlignment swfast(const string& seq1, const string& seq2, double al_m, double al_b, int width)
{
    int n0 = seq1.size();
    int n1 = seq2.size();
    
    // m and b are given such that i = m*j + b
    // first, calculate intercepts, so we know where to start and stop
    
    int j0 = (int)floor((-width/2 - al_b)/al_m);
    int j1 = (int)floor((n0+width/2 - al_b)/al_m);
    
    if (j0 < 0) j0 = 0;
    if (j0 >= n1) j0 = n1-1;
    if (j1 < 2) j1 = 2;
    if (j1 > n1) j1 = n1;
    
    // next, allocate arrays, with j0...j1 inclusive
    int* scores = new int[(j1-j0+1)*width]();
    uint8_t* steps = new uint8_t[(j1-j0+1)*width]();
    // and the starting i-indices into above arrays
    int* i0s = new int[(j1-j0+1)];
    // which we then calculate, even if out of range
    for (int j=j0; j<=j1; j++)
        i0s[j-j0] = (int)floor(al_m*j+al_b) - width/2;
    
    // what and where is the max score
    int maxScore = 0;
    int maxI = 0;
    int maxJ = 0;
    
    // go through all relevant columns and calculate
    // (skip the first column to make some if statements easier)
    for (int j=j0+1; j<=j1; j++)
    {
        // so j has range minimum 1 maximum n1 here
        
        // range for the current stripe
        int i0 = i0s[j-j0];
        int i1 = i0 + width - 1;
        // and the stripe to the left
        int p0 = i0s[j-j0-1];
        int p1 = p0 + width - 1;
        
        // clip the values of i0 to 1...n0
        if (i0 < 1) i0 = 1;
        if (i0 > n0) i0 = n0;
        if (i1 < 1) i1 = 1;
        if (i1 > n0) i1 = n0;
        
        // now let's go through and write each output elt.
        // first, create a fake pointer to this column's base
        // and previous column's base
        int* curscore = scores + (j-j0)*width - i0s[j-j0];
        int* prevscore = scores + (j-j0-1)*width - i0s[j-j0-1];
        uint8_t* curstep = steps + (j-j0)*width - i0s[j-j0];
        
        // now loop through i0,i1, which should totes work
        for (int i=i0; i<=i1; i++)
        {
            // best score
            int score=0;
            // best step
            uint8_t step=0;
            
            // note, never gonna have implicit insertions, since those
            // are always negative
            
            // check step from left
            if (i >= p0 && i <= p1)
            {
                int s = prevscore[i] + score_insert;
                if (s > score)
                {
                    score = s;
                    step = 1;
                }
            }
            // and above
            if (i > i0)
            {
                int s = curscore[i-1] + score_insert;
                if (s > score)
                {
                    score = s;
                    step = 2;
                }
            }
            // and above-left, match/mismatch
            if (i > p0 && i <= p1)
            {
                int s = prevscore[i-1] + 
                        ((seq1[i-1]==seq2[j-1])?score_match:score_mismatch);
                if (s >= score)
                {
                    score = s;
                    step = 3;
                }
            }
            else
            {
                // implicit match/mismatch, flag it as such
                int s = (seq1[i-1]==seq2[j-1])?score_match:score_mismatch;
                if (s >= score)
                {
                    score = s;
                    step = 255;
                }
            }
            
            curscore[i] = score;
            curstep[i] = step;
            
            if (score > maxScore)
            {
                maxScore = score;
                maxI = i;
                maxJ = j;
            }
        }
    }
    
    // now run backtrace
    int i = maxI;
    int j = maxJ;
    
    // create output struct
    SWAlignment align;
    align.score = maxScore;
    vector<int>& inds1 = align.inds1;
    vector<int>& inds2 = align.inds2;
    
    int nmatch = 0;
    
    while (i > 0 && j > 0)
    {
        int curscore = scores[(j-j0)*width - i0s[j-j0] + i];
        uint8_t curstep = steps[(j-j0)*width - i0s[j-j0] + i];
        
        // hit a 0?
        if (curscore <= 0)
            break;
                
        switch (curstep)
        {
            case 1:
                // step from left, j
                inds1.push_back(0);
                inds2.push_back(j);
                j--;
                break;
            case 2:
                // above, i
                inds1.push_back(i);
                inds2.push_back(0);
                i--;
                break;
            case 3:
                // les deux, match (or mismatch)
                inds1.push_back(i);
                inds2.push_back(j);
                if (i>0 && j>0 && seq1[i-1] == seq2[j-1])
                    nmatch++;
                i--;
                j--;
                break;
            case 255:
                // or just done, no previous step
                inds1.push_back(i);
                inds2.push_back(j);
                i = 0;
                j = 0;
                break;
            case 0:
            default:
                cerr << "urghhhh bad bad bad day" << endl;
        }
    }
    
    // reverse the vectors
    reverse(inds1.begin(),inds1.end());
    reverse(inds2.begin(),inds2.end());
    
    align.accuracy = 100.0*nmatch / (double)inds1.size();

    // and clean up
    delete scores;
    delete steps;
    delete i0s;
    
    return align;
}

SWAlignment swfull(const string& seq1, const string& seq2)
{
    int n1 = seq1.size();
    int n2 = seq2.size();
    
    // next, allocate full, initialized arrays
    int* scores = new int[(n1+1)*(n2+1)]();
    uint8_t* steps = new uint8_t[(n1+1)*(n2+1)]();
    
    // what and where is the max score
    int maxScore = 0;
    int maxI = 0;
    int maxJ = 0;
    
    // go through all relevant columns and calculate
    // (skip the first column to make some if statements easier)
    for (int j=1; j<=n2; j++)
    {
        // now let's go through and write each output elt.
        // first, create a fake pointer to this column's base
        // and previous column's base
        int* curscore = scores + j*(n1+1);
        int* prevscore = scores + (j-1)*(n1+1);
        uint8_t* curstep = steps + j*(n1+1);
        
        // now loop through the rows
        for (int i=1; i<=n1; i++)
        {
            // best score
            int score=0;
            // best step
            uint8_t step=0;
            int s;

            // check step from left
            s = prevscore[i] + score_insert;
            if (s > score)
            {
                score = s;
                step = 1;
            }
            s = curscore[i-1] + score_insert;
            if (s > score)
            {
                score = s;
                step = 2;
            }
            // and above-left, match/mismatch
            s = prevscore[i-1] + 
                    ((seq1[i-1]==seq2[j-1])?score_match:score_mismatch);
            if (s >= score)
            {
                score = s;
                step = 3;
            }
            
            curscore[i] = score;
            curstep[i] = step;
            
            if (score > maxScore)
            {
                maxScore = score;
                maxI = i;
                maxJ = j;
            }
        }
    }
    
    // now run backtrace
    int i = maxI;
    int j = maxJ;
    
    // create output struct
    SWAlignment align;
    align.score = maxScore;
    vector<int>& inds1 = align.inds1;
    vector<int>& inds2 = align.inds2;
    
    int nmatch = 0;
    
    while (i > 0 && j > 0)
    {
        int curscore = scores[j*(n1+1) + i];
        uint8_t curstep = steps[j*(n1+1) + i];
        
        // hit a 0?
        if (curscore <= 0)
            break;
                
        switch (curstep)
        {
            case 1:
                // step from left, j
                inds1.push_back(0);
                inds2.push_back(j);
                j--;
                break;
            case 2:
                // above, i
                inds1.push_back(i);
                inds2.push_back(0);
                i--;
                break;
            case 3:
                // les deux
                inds1.push_back(i);
                inds2.push_back(j);
                if (i>0 &&j >0 && seq1[i-1] == seq2[j-1])
                    nmatch++;
                i--;
                j--;
                break;
            case 0:
            default:
                cerr << "Smith-Waterman memory error" << endl;
        }
    }
    
    // reverse the vectors
    reverse(inds1.begin(),inds1.end());
    reverse(inds2.begin(),inds2.end());
    
    align.accuracy = 100.0*nmatch / (double)inds1.size();

    // and clean up
    delete scores;
    delete steps;
    
    return align;
}

SWAlignment fillinds(SWAlignment align)
{
    // struct passed by value so that we make a copy
    vector<int>& inds1 = align.inds1;
    vector<int>& inds2 = align.inds2;
    
    int i1 = inds1[0];
    int i2 = inds2[0];
    
    for (int i=0; i<inds1.size(); i++)
    {
        if (inds1[i] > 0)
            i1 = inds1[i];
        else
            inds1[i] = i1;
        
        if (inds2[i] > 0)
            i2 = inds2[i];
        else
            inds2[i] = i2;
    }
    
    return align;
}