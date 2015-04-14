/*==========================================================
 *
 * EventData.h
 *
 * Classes to hold event and model data and stuff.
 *
 *========================================================*/

#ifndef _EVENTDATA_H_
#define _EVENTDATA_H_

#include <vector>
#include <algorithm>

#include "AlignUtil.h"
#include "Sequence.h"

using namespace std;

// internal model struct
struct ModelData
{
    // supplied values
    double lev_mean[N_STATES];
    double lev_stdv[N_STATES];
    double sd_mean[N_STATES];
    double sd_stdv[N_STATES];

    // and additional skip-stay params
    double lik_skip;
    double lik_stay;
    double lik_extend;
    double lik_insert;
    
    // computed values
    double log_lev[N_STATES];
    double log_sd[N_STATES];
    // inverse gaussian parameters
    double sd_lambda[N_STATES];
    double log_lambda[N_STATES];
    
    // is it a complement model?
    bool complement;

    // set the model data from the outside
    // this is in a function so if we change what we have, we don't forget things
    // (because then it will stop compiling)
    void setData(double* levmean, double* levstdv, double* sdmean, double* sdstdv, bool complement)
    {
       // copy over and also precalculate logarithmed versions
        for (int i=0; i<N_STATES; i++)
        {
            lev_mean[i] = levmean[i];
            lev_stdv[i] = levstdv[i];
            sd_mean[i] = sdmean[i];
            sd_stdv[i] = sdstdv[i];
            log_lev[i] = log(levstdv[i]);
            log_sd[i] = log(sdstdv[i]);
        
            // oxford's scaling to get inverse gaussian parameters
            sd_lambda[i] = pow(sd_mean[i],3)/pow(sd_stdv[i],2);
            log_lambda[i] = log(sd_lambda[i]);
        }
    }

    // same comment as above
    void setParams(double prob_skip, double prob_stay, double prob_extend, double prob_insert)
    {
        lik_skip = log(prob_skip);
        lik_stay = log(prob_stay);
        lik_extend = log(prob_extend);
        lik_insert = log(prob_insert);
    }
};


// struct for an event
struct EventData
{
    // internal model structure corresponding to event
    ModelData model;
    // and the sequence data associated with the event
    // (the 2d or 1d original sequence)
    Sequence sequence;
    
    // how many levels?
    int length;
    // beginning and end of reference
    int refstart;
    int refend;
    
    // ref_align:   0 before first aligned position or after last aligned
    //             -1 at inserted levels between aligned bases
    //              or just index of aligned base otherwise
    
    // ref_like:    likelihood change from previous aligned state to here
    
    // ref_index:   like ref_align, except linearly interpolated instead
    //              of having 0 or -1; monotonic and used for fast indexing
    //              (and is generated from ref_align)

    vector<double> mean;
    vector<double> stdv;
    vector<double> log_stdv;
    vector<double> ref_align;
    vector<double> ref_index;
    vector<double> ref_like;
    
    // populate ref_index densely with indices and such from ref_align
    void updaterefs()
    {
        double al_m, al_b;  // alignment coefficients to extend past start
        int ra0, ra1;       // indices in ref_align
        
        ra0 = -1;
        ra1 = -1;
        refstart = -1;
        refend = -1;
        // ok, find first and last nonnegative index in ref_align
        for (ra0=0; ra0<this->length; ra0++)
            if (ref_align[ra0] > 0)
                break;
        for (ra1=this->length-1; ra1>=0; ra1--)
            if (ref_align[ra1] > 0)
                break;

        // if we don't have any, abort
        if (ra0 == this->length || ra1 < 0)
        {
            ref_index = vector<double>();
            return;
        }
        
        // now set ref start and end
        refstart = ref_align[ra0];
        refend = ref_align[ra1];
        
        // otherwise, initialize array as a copy
        ref_index = ref_align;
        
        // calculate linear alignment and stuff
        al_m = (ref_align[ra1]-ref_align[ra0])/(double)(ra1-ra0);
        al_b = ref_align[ra0] - al_m*ra0;
        // now fill in all of it
        int lastal = -1; // previous ref_align-populated state
        for (int i=0; i<this->length; i++)
        {
            // before first and last aligned states
            if (i < ra0 || i > ra1)
            {
                // just use linear extrapolation
                ref_index[i] = al_m*i + al_b;
            }
            // we have an alignment
            else if (ref_align[i] > 0)
            {
                if (lastal > 0)
                {
                    // loop from lastal + 1
                    double m = (ref_align[i]-ref_align[lastal])/(i-lastal);
                    // and fill in all the ref_indexes
                    for (int j=lastal+1; j<i; j++)
                        ref_index[j] = m*(j-lastal) + ref_align[lastal];
                }
                lastal = i;
                // ref_index[i] already is ref_align[i]
            }
        }
    }
    
    // find just the one single index
    int getrefstate(int refind)
    {
        // if empty, just play dumb
        if (ref_index.size() == 0)
            return 0;
        // find first position greater than or equal to refind
        vector<double>::iterator it = std::lower_bound(ref_index.begin(),ref_index.end(),refind);
        // if refind is less than ref_index[0], this will be fine, since it == ref_index.begin() then
        // and if refind is more than the end of ref_index, also fine, since default behavior is to return
        // ref_index.end() if not found
        return it - ref_index.begin();
    }

    // find all the indices that lie within a certain ref index
    // (only if it exists, yadda yadda)
    vector<int> getrefstates(int refind)
    {
        vector<int> inds;
        
        // find index of refind in ref_index array
        vector<double>::iterator it = std::find(ref_index.begin(),ref_index.end(),refind);
        // if we didn't find it, return nada
        if (it == ref_index.end())
            return inds;
        
        int i = it - ref_index.begin();
        inds.push_back(i);
        for (i++; i<length && ref_align[i] <= refind; i++)
            if (ref_align[i] > 0)
                inds.push_back(i);
        
        return inds;
    }

    // this function works out well for both Matlab and Cython
    // (note: cython is not great with non-default constructors, so...)
    void setData(int count, double* meanpr, double* stdvpr, double* ref_alignpr, double* ref_likepr)
    {
        this->length = count;
        
        mean = vector<double>(meanpr,meanpr+length);
        stdv = vector<double>(stdvpr,stdvpr+length);
        ref_align = vector<double>(ref_alignpr,ref_alignpr+length);
        ref_like = vector<double>(ref_likepr,ref_likepr+length);
        
        // populate the log-stdv vector for precomputation
        log_stdv = vector<double>(this->length);
        for (int i=0; i<this->length; i++)
            this->log_stdv[i] = log(this->stdv[i]);
        
        // then initialize the internal indexing helpers
        updaterefs();
    }

    
    EventData() : length(0)
    {}
};


#endif
