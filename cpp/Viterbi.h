/*==========================================================
 *
 * Viterbi.h
 *
 * Helper functions/transition matrix for 1d, 2d C++ Viterbi algorithms
 * The Viterbi code is older code, barely used, and pretty junky...
 *
 *========================================================*/

#ifndef _VITERBI_H_
#define _VITERBI_H_

#include <stdint.h>
#include <vector>

#include "EventData.h"
#include "AlignUtil.h"
#include "Sequence.h"

using namespace std;

// previous state macro/inline function
inline int prev_state(int state, int ind)
{ return (state >> 2)+(ind<<8); }
inline int next_state(int state, int ind)
{ return ((state << 2)&(N_STATES-1))+ind; }
// and with multiple advances
inline int prev_state(int state, int ind, int nsteps)
{ return (state >> (2*nsteps))+(ind<<(10-2*nsteps)); }
inline int next_state(int state, int ind, int nsteps)
{ return ((state << (2*nsteps))&(N_STATES-1))+ind; }

// return base at ind, ind = 0...4, 0 is leftmost, 4 is rightmost
inline char get_base(int state, int ind)
{
    static char* bases = "ACGT";
    return bases[3 & (state>>(2*(4-ind)))];
}

// returns rev. complement of a 5-mer index
inline int complement_state(int state)
{
    int comp = 0;
    
    // reverse two bits at a time, taking xor of said bits
    for (int i=0; i<5; i++)
    {
        comp <<= 2;
        comp += (state&3)^3;
        state >>= 2;
    }
    return comp;
}

// normalize an N_STATES vector
inline void normvec(double* vec)
{
    double tot = 0;
    for (int i=0; i<N_STATES; i++)
        tot += vec[i];
    tot = 1.0/tot;
    for (int i=0; i<N_STATES; i++)
        vec[i] *= tot;
}

// the actual function to generate mutations
vector<Sequence> ViterbiMutate(vector<EventData>& events, int nkeep, 
        double skip_prob, double stay_prob, double mut_min, double mut_max, bool verbose);

// a helper function to convert a vector of states into a sequence class
Sequence StatesToSequence(vector<int> states);

#endif
