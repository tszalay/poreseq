/*==========================================================
 *
 * Alignment.h
 *
 * Class for alignment of a single event/model pair with a given sequence.
 *
 *========================================================*/

#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include "AlignUtil.h"
#include "EventData.h"
#include "Sequence.h"

#include <stdint.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

using namespace std;


// how many likelihood matrices there are
const int AL_LIKS=2;

/**************** NOTE NOTE NOTE NOTE NOTE ON INDEXING!!!!!
 * i and j (and all variants thereof) refer to the actual index in an
 * implicit matrix of size (n0+1)*(n1+1), where the first row and
 * column are all 0.
 */


// this is a struct for holding data along with a score
// and keeping the data corresponding to best score
struct MaxInfo
{
    double score;
    int i;
    int j;
    
    MaxInfo() : score(0), i(0), j(0)
    {}
    
    MaxInfo(double _s, int _i, int _j) : score(_s), i(_i), j(_j)
    {}
};

// assign max of the two to LHS, and track i,j
MaxInfo& operator<<= (MaxInfo& mi1, const MaxInfo& mi2);

// struct for saving a partial column of a sparse matrix only
// (this is designed to be a shared ptr, so no copy)
struct AlignColumn
{
    const int i0;             // starting row of the data
    const int col;            // column of the data
    const int length;         // how many data points saved?
    double* liks;             // pointer to the scores array, M x NUM
    double* probs;            // observation probabilities
    uint8_t* steps;           // and the backpointers
    MaxInfo maxScore;         // best score/i/j in matrix up to this point
    
    // the actual data, allocated as one chunk, to save on malloc time
    double* _data;
    
    // create an empty one
    AlignColumn(int num, int ind0, int ind1) : i0(ind0), col(ind1), length(num)
    {
        // parens at the end make this equivalent to calloc
        this->_data = new double[this->length*(AL_LIKS+1)+(this->length*AL_LIKS)/8+1]();
        this->liks = this->_data;
        this->probs = this->_data + this->length*AL_LIKS;
        this->steps = (uint8_t*)(this->_data + this->length*(AL_LIKS+1));
    }
    
    // get a pointer to an element in likelihood array
    // the "arr" index is for which array, of which there are AL_LIKS
    // (in practice, two, one for normal and one for stays)
    double* getPointer(int ind, int arr)
    {
        return this->liks + arr*this->length + (ind-this->i0);
    }
    
    // get a pointer to an element in backpointer array
    uint8_t* getStep(int ind, int arr)
    {
        return this->steps + arr*this->length + (ind-this->i0);
    }

    // get a pointer to an element in observation array
    double* getProb(int ind)
    {
        return this->probs + (ind-this->i0);
    }

    
    ~AlignColumn()
    {
        delete _data;
    }
};


typedef shared_ptr<AlignColumn> AlignPointer;


class Alignment
{
    // references to external data
    EventData* event;
    Sequence* sequence;
    AlignParams params;
    
    // our internal vectors containing the columns
    // both from left-to-right and right-to-left
    vector<AlignPointer> scores;
    vector<AlignPointer> scores_back;

    // what stripe width we're currently using
    int stripe_width;
    
public:
    
    double getMax()
    {
        return max(scores.back()->maxScore.score,scores_back.back()->maxScore.score);
    }
    
    // construct with an event, model, sequence, etc.
    Alignment(Sequence& seq, EventData& ed, AlignParams& par);
    
    // refill all the columns from first empty one, for a new sequence
    // (also runs a backtrace)
    void update(Sequence& seq);
    
    // clear all columns, in proper order
    void clear();
    
    // fill in the rest of the columns, forwards
    void fillColumns();
    
    // fill in the rest of the columns, backwards
    void fillColumnsBack();
    
    // fill in a few coluns
    void fillColumns(int n);
    
    // fill in the next column of the alignment, if more remain, and
    // return the max value recorded in that column
    double fillColumn();
    
    // fill in the next column of the alignment, if more remain, and
    // return the max value recorded in that column
    double fillColumnBack();
    
    // calculate (not estimate) the score from a given mutation
    double scoreMutation(const MutInfo& mut, Sequence& mutseq);
    
    // run a backtrace and populate an event's ref_align, if given
    EventData* backtrace();
    


private:
    
    double columnMax(int refind)
    {
        // calculate single column maximum
        int raf = refind;
        int rab = sequence->states.size()-refind+1;
        
        return columnMax(raf, rab);
    }
    
    // calculate column max score, with fwd and back indices
    // this is the function that uses the combination of a column from
    // each of the forward and back matrices to quickly get the likelihood
    double columnMax(int raf, int rab)
    {
        // basic bounds checking, just cuz, this can happen
        if (raf >= scores.size()) raf = scores.size()-1;
        if (rab >= scores_back.size()) rab = scores_back.size()-1;
        if (raf < 0) raf = 0;
        if (rab < 0) rab = 0;

        // calculate max product of forward and back in each column
        double sm = 0;
        for (int jf=1; jf<=event->length; jf++)
        {
            // pointers to columns in fwd and back
            AlignPointer& sf = scores[raf];
            AlignPointer& sb = scores_back[rab];
            int jb = event->length-jf+1;

            // now add together
            for (int k=0; k<2; k++)
            {
                double s = 0;
                // make sure we check ranges
                if (jf >= sf->i0 && jf < sf->i0+sf->length)
                    s += *sf->getPointer(jf,k);
                if (jb >= sb->i0 && jb < sb->i0+sb->length)
                    s += *sb->getPointer(jb,k);

                sm = max(s,sm);
            }
            sm = max(sm, sf->maxScore.score);
            sm = max(sm, sb->maxScore.score);
        }
        return sm;
    }
};

#endif