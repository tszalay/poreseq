/*==========================================================
 *
 * Viterbi.cpp
 *
 * Contains Viterbi mutation functions
 *
 *========================================================*/

#include <stdint.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

#include "Viterbi.h"


using namespace std;


// holds viterbi likelihoods etc etc
struct V_LIK
{
    double liks[N_STATES];              // accumulated likelihood so far, best-path
    int backptrs[N_STATES];             // viterbi backpointers
    double fwdprobs[N_STATES];          // forward probability values
    
    // construct from previous likelihood struct
    V_LIK(unique_ptr<V_LIK>& prevlik, vector<double>& obs, double skip_prob, double stay_prob);
    V_LIK() {};
    
    // get a random backpointer, for a given state, using a given
    // transition matrix (not log-form)
    int randbp(int curstate, double atten, const vector<double>& T);
    
};
typedef unique_ptr<V_LIK> V_PTR;

V_LIK::V_LIK(V_PTR& prevlik, vector<double>& obs, double skip_prob, double stay_prob)
{
    // number of skips to consider ('1 skip' is a normal step)
    const int nskip = 3;
    const double skip_lik = log(skip_prob);
    const double stay_lik = log(stay_prob);
    
    // all destination states
    for (int curst=0; curst<N_STATES; curst++)
    {
        // best numbers for this state
        double maxlik = -inf;
        int maxptr = -1;
        // forward algorithm probability
        double fwdprob = 0.0;

        // a loop for each number of skips
        double sp = 0.25;
        double lsp = log(0.25);
        
        for (int j=1; j<=nskip; j++)
        {
            // loop through all possible states we came from
            for (int k=0; k < 1<<(2*j); k++)
            {
                int prevst = prev_state(curst, k, j);
                
                double l = obs[curst] + lsp;
                l += prevlik->liks[prevst];
                fwdprob += sp * prevlik->fwdprobs[prevst];
                if (l > maxlik)
                {
                    maxlik = l;
                    maxptr = prevst;
                }
            }
            // and change prob for next time through
            sp = sp*0.25*skip_prob;
            lsp = lsp + log(0.25) + skip_lik;
        }
        
        // now finally check stays, add it to self-transition prob
        double l = obs[curst]+stay_lik+prevlik->liks[curst];
        if (l > maxlik)
        {
            maxlik = l;
            maxptr = curst;
        }

        fwdprob += stay_prob * prevlik->fwdprobs[curst];

        // multiply forward probability by the observation probability
        fwdprob *= exp(obs[curst]);

        // and partially sort the array to get best few likelihoods
        // save the best likelihood only into lik array
        this->liks[curst] = maxlik;
        this->backptrs[curst] = maxptr;
        this->fwdprobs[curst] = fwdprob;
    }
    
    // normalize forward probabilities
    normvec(this->fwdprobs);
}

// generate a random backpointer using forward probs (makes sense, I think?)
int V_LIK::randbp(int curstate, double atten, const vector<double>& T)
{
    // random value 0...1
    double r = rand() / (double(RAND_MAX)+1);
    // and value to compare to
    double cumsum = 0;
    // now renormalize probs multiplied by transitions
    double probs[N_STATES];
    
    for (int i=0; i<N_STATES; i++)
        probs[i] = T[i + curstate*N_STATES] * pow(this->fwdprobs[i],atten);

    // and normalize it
    normvec(probs);
    
    // and find the one that is hit by r
    for (int i=0; i<N_STATES; i++)
    {
        cumsum += probs[i];
        if (r < cumsum)
            return i;
    }
    // if we didn't find one yet, just some numerical slop, it should
    // be the last state
    // (this doesn't happen in practice)
    return N_STATES-1;
}

// create a transition matrix
vector<double> buildT(double skip_prob, double stay_prob)
{
    // initialize to 0
    vector<double> T(N_STATES*N_STATES,0);
        
    // number of skips to consider ('1 skip' is a normal step)
    const int nskip = 4;
    
    // all destination states
    for (int curst=0; curst<N_STATES; curst++)
    {
        // current row (or column, or whatever)
        // small steps in T represent origin states
        auto Tcol = T.begin() + N_STATES*curst;
        
        // a loop for each number of skips
        double sp = 0.25;
        for (int j=1; j<=nskip; j++)
        {
            // loop through all possible new states
            for (int k=0; k < 1<<(2*j); k++)
            {
                int prevst = prev_state(curst, k, j);
                Tcol[prevst] += sp;
            }
            // and change prob for next time through
            sp = sp*0.25*skip_prob;
        }
    }
    
     // now make sure transition matrix knows about stays properly
    for (int i=0; i<N_STATES; i++)
        T[i*(1+N_STATES)] = stay_prob;
    
    return T;
}

Sequence StatesToSequence(vector<int> states)
{
    string seq;
    
    int curstate = states[0];
    seq.push_back(get_base(curstate,0));
    
    int nmismatch = 0;
    int nskip = 0;
    int nstay = 0;
    
    // now loop through the rest of the states
    for (int i=1; i<states.size(); i++)
    {
        if (curstate == states[i])
        {
            // we stayed where we were, no transition
            // curstate stays, states stays, nothing added to sequence
            nstay = nstay + 1;
            continue;
        }
        // otherwise, we have a new state, figure out what the possible
        // transition could have been
        // (nskips == 1 is a regular step)
        for (int nskips=1; nskips<=4; nskips++)
        {
            for (int ind=0; ind < (1<<(2*nskips)); ind++)
            {
                if (next_state(curstate,ind,nskips)==states[i])
                {
                    // we found it! first, insert missing letters
                    // already added curstate's zeroth base
                    for (int j=1; j<=nskips; j++)
                        seq.push_back(get_base(curstate,j));
                    
                    // now update curstate
                    curstate = states[i];
                    // and write relevant vars, '1' skip is 0
                    // i really apologize for this line of code
                    // the naming and indexing both suck
                    nskip = nskip + nskips - 1;
                    break;
                }
            }
            // we found it in inner loop
            if (curstate==states[i])
                break;
            
        }
        // didn't find it?
        if (curstate != states[i])
        {
            // likely we have a mismatch
            // just do a sudden transition
            curstate = states[i];
            // take the middle letter of it
            seq.push_back(get_base(curstate,0));
            nmismatch = nmismatch+1;
        }
    }
    
    // now pop on the last bases we didn't do yet
    for (int i=1; i<=4; i++)
        seq.push_back(get_base(curstate,i));
    
    return Sequence(seq);
}

vector<Sequence> ViterbiMutate(vector<EventData>& events, int nkeep, 
        double skip_prob, double stay_prob, double mut_min, double mut_max, bool verbose)
{
    // contain the dynamically growing array of scores
    vector<V_PTR> scores;

    int n_seq = events.size();

    // let's make a blank, dummy initial likelihood
    scores.push_back(V_PTR(new V_LIK()));
    V_PTR& lik0 = scores.back();
        
    for (int i=0; i<N_STATES; i++)
    {
        lik0->liks[i] = 0;
        lik0->backptrs[i] = -1;
        lik0->fwdprobs[i] = 1.0/N_STATES;
    }
	
    if (verbose)
        cerr << "Viterbi";
    
    // find event that has the lowest ref_align, and start there
    int refind = events[0].refstart;
    for (int i=0; i<n_seq; i++)
        refind = min(refind,events[i].refstart);

    vector<double> obs(N_STATES*n_seq);

    // now loop through each level
    while (1)
    {
        // calculate log-normalized observations
        std::fill(obs.begin(),obs.end(),0);
        
        // figure out which strands actually matter here
        int nlik = 0;
        for (int k=0; k<n_seq; k++)
        {
            vector<int> inds = events[k].getrefstates(refind);
            
            // does this strand have any states that line up with this ref state?
            if (inds.size() == 0)
                continue;
            
            nlik++;
            
            // now, average all these thingies, and stuff
            double lvl = 0;
            double sd = 0;
            for (int j=0; j<inds.size(); j++)
            {
                lvl += events[k].mean[inds[j]];
                sd += events[k].stdv[inds[j]];
            }
            lvl = lvl/inds.size();
            sd = sd/inds.size();
            
            ModelData& model = events[k].model;

            // calculate observation likelihoods for each state
            for (int j=0; j<N_STATES; j++)
            {
                double l = lognormpdf(lvl,model.lev_mean[j],model.lev_stdv[j],model.log_lev[j]);
                l += logigpdf(sd,model.sd_mean[j],model.sd_lambda[j],
                                log(sd),model.log_lambda[j]);
                obs[j*n_seq + nlik-1] = l;
            }
        }
        
        // how many align here?
        int nalhere = 0;
        for (int k=0; k<n_seq; k++)
            if (refind >= events[k].refstart && refind <= events[k].refend)
                nalhere++;
        
        // if too few aligned strands, skip over this state (or end, if done)
        if (nlik <= nalhere*0.2)
        {
            // stop, if no more strands
            if (nalhere == 0)
                break;
            
            // or just skip this one
            refind++;
            continue;
        }
        
        if (nlik > 1)
        {
            // sort obs likelihoods, ascendingally
            for (int j=0; j<N_STATES; j++)
                std::sort(obs.begin()+j*n_seq,obs.begin()+j*n_seq+nlik);
            // and then calculate means likelihood of the last few elements only
            int nskip = floor(nlik*0.25);
            if (nskip > nlik-2) nskip = 0;
            for (int j=0; j<N_STATES; j++)
            {
                double lik = 0.0;
                for (int k=nskip; k<nlik; k++)
                    lik += obs[j*n_seq + k];
                // now put it in the first column, this is shitty code
                obs[j] = lik / (nlik - nskip);
            }
        }
        else
        {
            // nothing really to do
            for (int j=0; j<N_STATES; j++)
                obs[j] = obs[j*n_seq];
        }

        // save a new likelihood
        scores.push_back(V_PTR(new V_LIK(scores.back(),obs,skip_prob,stay_prob)));        

        // and advance
        refind++;
        
        if (refind%200 == 0)
        {
            if (verbose)
            {
                cerr << ".";
                cerr.flush();
            }
        }
    }

    if (verbose)
    {
        cerr << "\n";
        cerr.flush();
    }

    // outputs
    vector<Sequence> seqs;
    vector<int> states;

    // run backtrace on main path
    // find largest final likelihood
    double* mlik = std::max_element(scores.back()->liks,scores.back()->liks+N_STATES);
    int startst = mlik - scores.back()->liks;
    int curst = startst;
    
    int n = scores.size() - 1;
    
    if (nkeep == 0)
    {
        for (int i=n-1; i>=0; i--)
        {
            states.push_back(curst);
            curst = scores[i+1]->backptrs[curst];
        }
        
        // they were backwards, so flip them
        std::reverse(states.begin(),states.end());
        
        // and turn them into a sequence
        seqs.push_back(StatesToSequence(states));

        return seqs;
    }

    // now for the individual ones, create transition matrix for back-stepping
    vector<double> T = buildT(skip_prob, stay_prob);
    
    for (int k=0; k<nkeep; k++)
    {
        // put all the states traversed into states variable
        states.clear();
        curst = startst;
        
        // now backtrack
        for (int i=n-1; i>=0; i--)
        {
            states.push_back(curst);
            curst = scores[i+1]->randbp(curst,mut_min+(mut_max-mut_min)*k/(double)nkeep,T);
        }
        
        // they were backwards, so flip them
        std::reverse(states.begin(),states.end());
        
        // and turn them into a sequence
        seqs.push_back(StatesToSequence(states));
    }
    
    return seqs;
}
