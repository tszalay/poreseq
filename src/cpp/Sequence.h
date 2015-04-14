/*==========================================================
 *
 * Sequence.h
 *
 * Class for manipulating/mutating/whatever a sequence.
 * This implementation is slow, will change if necessary.
 *
 *========================================================*/

#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include "AlignUtil.h"
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

struct Sequence
{
    string          bases;  // the ACGT sequence
    vector<int>     states; // the states corresponding to above
    
    // note that this means states is 4 elements shorter than bases
    // which is kind of necessary, but a bit unfortunate
    
    Sequence() {}
    
    // construct it from a C++ string
    Sequence(string seq) : bases(seq)
    {
        this->populateStates();
    }
    
    // construct using a specified mutation
    Sequence(const Sequence& original, const MutInfo& mut)
    {
        // if mutation is past end, copy-construct 
        if (mut.start >= original.bases.size())
        {
            this->bases = original.bases;
            this->states = original.states;
            return;
        }

        // first, copy over the original states
        this->bases = original.bases.substr(0,mut.start);
        // then insert the new ones
        this->bases += mut.mut;
        // and then the remaining ones, to end of string
        int remind = mut.start+mut.orig.size();
        if (remind < original.bases.size())
            this->bases += original.bases.substr(remind);

        // then recompute the state indices
        this->populateStates();
    }

private:

    // this function takes the string of bases and turns it into
    // an array of states 0...1023
    void populateStates()
    {
        // only do this if we have enough bases
        if (bases.size() < 5)
            return;
        
        // replace letter bases with integers for next step
        string baseint = bases;
        std::replace(baseint.begin(), baseint.end(), 'A', char(0));
        std::replace(baseint.begin(), baseint.end(), 'C', char(1));
        std::replace(baseint.begin(), baseint.end(), 'G', char(2));
        std::replace(baseint.begin(), baseint.end(), 'T', char(3));

        // initialize it with the first few bases
        int curstate = 0;
        for (int i=0; i<4; i++)
            curstate = (curstate << 2) + baseint[i];
        
        // and then add in the rest
        for (int i=4; i<baseint.size(); i++)
        {
            // do we have an invalid state influence?
            if (baseint[i-4] < 4)
            {
                // no, proceed as normal
                curstate = (N_STATES-1)&((curstate << 2) + baseint[i]);
                this->states.push_back(curstate);
            }
            else
            {
                // yes, there is a '-' that would influence this state
                curstate = 0;
                this->states.push_back(-1);
            }
        }
    }
};

#endif