/*==========================================================
 *
 * Alignment.cpp
 *
 * Contains the "guts" of the alignment calculation, this is where it gets
 * a bit on the technical side. Not all of the code here is used, and
 * there is a lot of duplication in the forward and backward computation,
 * so tread lightly.
 *
 *========================================================*/

#include "Alignment.h"

using namespace std;

// how many backpointer possibilities we have
const int AL_BACKS=6;
// and their values/order, for switches
enum AL_MOVE
{
    L_SKIP = 0,
    UL_MATCH = 1,
    U_INSERT = 2,
    UL_IGNORE = 3,
    U_STAY = 4,
    U_EXTEND = 5,
    Z_IMPLICIT = 255
};

// assign max of the two to LHS, and track i,j
MaxInfo& operator<<= (MaxInfo& mi1, const MaxInfo& mi2)
{
    if (mi2.score > mi1.score)
        mi1 = mi2;
    return mi1;
}

Alignment::Alignment(Sequence& seq, EventData& ed, AlignParams& par) : event(&ed),
            sequence(&seq), params(par)
{
    // create a blank row
    scores.push_back(AlignPointer(new AlignColumn(event->length+1,0,0)));
    scores_back.push_back(AlignPointer(new AlignColumn(event->length+1,0,0)));

    // and set the default stripe width to realign
    this->stripe_width = params.realign_width;
    
    // the blank row will have a default MaxInfo that is suitable

    // if there is no existing reference alignment, don't use this event
    if (event->ref_index.size() == 0)
    {
        if (params.verbose)
        {
            cerr << "F";
            cerr.flush();
        }
        this->stripe_width = 0;
    }
}

// clear and re-populate sparse matrix with new sequence
void Alignment::update(Sequence& seq)
{
    this->clear();
    
    this->sequence = &seq;
    
    // fill forward and backward
    fillColumns();
    fillColumnsBack();
    backtrace();
}

void Alignment::clear()
{
    while (scores.size() > 1)
        scores.pop_back();
    while (scores_back.size() > 1)
        scores_back.pop_back();
}

void Alignment::fillColumns()
{
	if (this->stripe_width == 0)
		return;

    while (scores.back()->col < sequence->states.size())
        fillColumn();
}

void Alignment::fillColumnsBack()
{
	if (this->stripe_width == 0)
		return;

    while (sequence->states.size()+scores_back.back()->col > 0)
        fillColumnBack();
}

// fill the next n columns (if we need to)
void Alignment::fillColumns(int n)
{
    for (int i=0; i<n; i++)
        fillColumn();
}

// fills the column immediately after the backmost column
// (if we can fill another one)
// returns a meaningless value that I don't think is used anywhere
double Alignment::fillColumn()
{
    // go on to the next one
    int refind = scores.back()->col + 1;
    if (refind > sequence->states.size())
        return 0;
		
	if (this->stripe_width == 0)
		return 0;


    // shorthand
    int n0 = event->length;
    
    int refstate = sequence->states[refind-1];
    
    int curwid = this->stripe_width;
    
    int imid = 1;
    // get index of closest ref_index
    if (event->ref_index.size() > 0)
        imid = event->getrefstate(refind);
    
    // if we are really searching before or after expected aligned region
    // (and not searching full-frame) let's speed it up a bit/lot
    if (curwid < n0 && (imid < -10 || imid > n0 + 10))
        curwid = 5;

    // check all limits and whatnot
    if (imid < 1) imid = 1;
    if (imid > n0) imid = n0;

    // upper and lower index to save, inclusive
    int i0 = imid-curwid;
    int i1 = imid+curwid;

    if (i0 < 1) i0 = 1;
    if (i1 > n0) i1 = n0;
    
    // now we know how much we're gonna fill in, allocate corresponding
    // scores vector and pull out their pointers
    scores.push_back(AlignPointer(new AlignColumn(i1-i0+1,i0,refind)));
    AlignPointer& curscore = scores[scores.size()-1];
    AlignPointer& prevscore = scores[scores.size()-2];
    
    
    // give curscore the best previous MaxInfo
    curscore->maxScore = prevscore->maxScore;
	
	// if we have an invalid state (eg. the sequence contains '-')
	// we can just return here, and save us the computation time
	if (refstate < 0)
		return curscore->maxScore.score;

    // write all the relevant observation likelihoods
    double* lik_obs = curscore->getProb(0);
    for (int i=i0; i<=i1; i++)
    {
        lik_obs[i] = lognormpdf(event->mean[i-1],event->model.lev_mean[refstate],event->model.lev_stdv[refstate],event->model.log_lev[refstate]);
        //lik_obs[i] += lognormpdf(event->stdv[i-1],event->model.sd_mean[refstate],event->model.sd_stdv[refstate],event->model.log_sd[refstate]);
        lik_obs[i] += logigpdf(event->stdv[i-1],event->model.sd_mean[refstate],event->model.sd_lambda[refstate],
                                event->log_stdv[n0-i],event->model.log_lambda[refstate]);
        lik_obs[i] += params.lik_offset;
    }
    
    // this is sketchy, because it points to nonexistant memory
    // but, in my defense, i'm, like, _really_ lazy
    // (so eg. curlik[0] is a segfault, since due to the sparsity only
    // curlik[i0]...curlik[i1] exist)
    double* curlik = curscore->getPointer(0,0);
    uint8_t* curstep = curscore->getStep(0,0);
    double* prevlik = prevscore->getPointer(0,0);
    // and pointers to gap/extension matrix, called "stay" matrix
    double* curstay = curscore->getPointer(0,1);
    uint8_t* curstaystep = curscore->getStep(0,1);
    double* prevstay = prevscore->getPointer(0,1);
    // initialize top element of stay matrix to -inf, not 0
    // since we can't start on a stay
    curstay[i0] = -inf;

    int p0 = prevscore->i0;
    int p1 = prevscore->i0+prevscore->length-1;

    for (int i=i0; i<=i1; i++)
    {
        double liks[AL_BACKS] = {0.0,0.0,0.0,0.0,-inf,-inf};
        uint8_t likbp[AL_BACKS] = {0,1,2,3,4,5};
        double lobs = lik_obs[i];

        // do implicit zeros if our indices don't work out
        if (i >= p0 && i <= p1)
        {
            // valid
            liks[L_SKIP] = prevlik[i] + event->model.lik_skip;
        }
        else
        {
            // implicit
            liks[L_SKIP] = event->model.lik_skip;
            likbp[L_SKIP] = Z_IMPLICIT;
        }

        if (i > p0 && i <= p1)
        {
            // valid
            liks[UL_MATCH] = prevlik[i-1] + lobs;
            // or ignore
            liks[UL_IGNORE] = prevlik[i-1] + event->model.lik_insert;
        }
        else
        {
            // implicit
            liks[UL_MATCH] = lobs;
            likbp[UL_MATCH] = Z_IMPLICIT;
        }

        // now here is where our stay matrices come into play
        // (no implicit transitions from above; can't stay or extend)
        if (i > i0)
        {
            // start a new stay
            liks[U_STAY] = curlik[i-1] + lobs + event->model.lik_stay;
            // or a new stay via insertion
            liks[U_INSERT] = curlik[i-1] + event->model.lik_insert;
            // or extend a stay, from stay matrix to stay matrix
            liks[U_EXTEND] = curstay[i-1] + lobs + event->model.lik_extend;
        }

        // first, update stay matrix with stays or extends
        for (int k=4; k<6; k++)
        {
            if (liks[k] > curstay[i])
            {
                curstay[i] = liks[k];
                curstaystep[i] = k;
            }
        }
        // note: stay matrix doesn't go into maxScore
        // just means we don't end on a stay - oh well...

        // ok, now first check transitions from main matrix
        for (int k=0; k<4; k++)
        {
            if (liks[k] > curlik[i])
            {
                curlik[i] = liks[k];
                // gotta do this to account for implicit zeros
                curstep[i] = likbp[k];
            }
        }

        // now check "transition" back from stay matrix
        if (curstay[i] > curlik[i])
        {
            curlik[i] = curstay[i];
            curstep[i] = U_STAY;
        }

        // and save, if it's the best
        curscore->maxScore <<= MaxInfo(curlik[i],i,refind);
    }
    
    return curscore->maxScore.score;
}

// same as previous function, but for backwards matrices
// so some things are reversed
// the #### signs more or less denote things that changed vs. the above fun
// except for the things that I probably forgot
double Alignment::fillColumnBack()
{
    // go on to the next one
    //#############################
    int colind = scores_back.back()->col - 1;
    int refind = sequence->states.size()+colind+1;
    if (refind <= 0)
        return 0;
		
	if (this->stripe_width == 0)
		return 0;

    
    // shorthand
    int n0 = event->length;
    
    int refstate = sequence->states[refind-1];
    
    int curwid = this->stripe_width;
    
    int imid = 1;
    // get index of closest ref_index (reversed)
    //#############################
    if (event->ref_index.size() > 0)
        imid = n0-event->getrefstate(refind)+1;
    
    
        
    // if we are really searching before or after expected aligned region
    // (and not searching full-frame) let's speed it up a bit/lot
    if (curwid < n0 && (imid < -10 || imid > n0 + 10))
        curwid = 5;

    // check all limits and whatnot
    if (imid < 1) imid = 1;
    if (imid > n0) imid = n0;

    // upper and lower index to save, inclusive
    int i0 = imid-curwid;
    int i1 = imid+curwid;

    if (i0 < 1) i0 = 1;
    if (i1 > n0) i1 = n0;
    
    // now we know how much we're gonna fill in, allocate corresponding
    // scores vector and pull out their pointers
    scores_back.push_back(AlignPointer(new AlignColumn(i1-i0+1,i0,colind)));
    //#############################
    AlignPointer& curscore = scores_back[scores_back.size()-1];
    AlignPointer& prevscore = scores_back[scores_back.size()-2];
    
    // give curscore the best previous MaxInfo
    curscore->maxScore = prevscore->maxScore;
	
    // if we have an invalid state (eg. the sequence contains '-')
	// we can just return here, and save us the computation time
	if (refstate < 0)
		return curscore->maxScore.score;

    // write all the relevant observation likelihoods
    double* lik_obs = curscore->getProb(0);
    for (int i=i0; i<=i1; i++)
    {
        //#############################
        // flip event too
        lik_obs[i] = lognormpdf(event->mean[n0-i],event->model.lev_mean[refstate],event->model.lev_stdv[refstate],event->model.log_lev[refstate]);
        //lik_obs[i] += lognormpdf(event->stdv[n0-i],event->model.sd_mean[refstate],event->model.sd_stdv[refstate],event->model.log_sd[refstate]);
        lik_obs[i] += logigpdf(event->stdv[n0-i],event->model.sd_mean[refstate],event->model.sd_lambda[refstate],
                                event->log_stdv[n0-i],event->model.log_lambda[refstate]);
        lik_obs[i] += params.lik_offset;
    }
    
    // this is sketchy, because it points to nonexistant memory
    // but, in my defense, i'm, like, _really_ lazy
    double* curlik = curscore->getPointer(0,0);
    uint8_t* curstep = curscore->getStep(0,0);
    double* prevlik = prevscore->getPointer(0,0);
    double* curobs = curscore->getProb(0);
    double* prevobs = prevscore->getProb(0);
    // and pointers to gap/extension matrix, called "stay" matrix
    double* curstay = curscore->getPointer(0,1);
    uint8_t* curstaystep = curscore->getStep(0,1);
    double* prevstay = prevscore->getPointer(0,1);
    // initialize top element of stay matrix to -inf, not 0
    // since we really super duper can't be there
    curstay[i0] = -inf;

    int p0 = prevscore->i0;
    int p1 = prevscore->i0+prevscore->length-1;

    for (int i=i0; i<=i1; i++)
    {
        double liks[AL_BACKS] = {0.0,0.0,0.0,0.0,-inf,-inf};
        uint8_t likbp[AL_BACKS] = {0,1,2,3,4,5};
        
        // do implicit zeros if our indices don't work out
        if (i >= p0 && i <= p1)
        {
            liks[L_SKIP] = prevlik[i] + event->model.lik_skip;
        }
        else
        {
            liks[L_SKIP] = event->model.lik_skip;
            likbp[L_SKIP] = Z_IMPLICIT;
        }

        if (i > p0 && i <= p1)
        {
            liks[UL_MATCH] = prevlik[i-1] + prevobs[i-1];
            liks[UL_IGNORE] = prevlik[i-1] + event->model.lik_insert;
        }
        else
        {
            liks[UL_MATCH] = 0;
            likbp[UL_MATCH] = Z_IMPLICIT;
        }

        // now here is where our stay matrices come into play
        // (no implicit transitions from above; can't stay or extend)
        if (i > i0)
        {
            // start a new stay
            liks[U_STAY] = curlik[i-1] + curobs[i-1] + event->model.lik_stay;
            // or a new stay via insertion
            liks[U_INSERT] = curlik[i-1] + event->model.lik_insert;
            // or extend a stay, from stay matrix to stay matrix
            liks[U_EXTEND] = curstay[i-1] + curobs[i-1] + event->model.lik_extend;
        }

        // first, update stay matrix with stays or extends
        for (int k=4; k<6; k++)
        {
            if (liks[k] > curstay[i])
            {
                curstay[i] = liks[k];
                curstaystep[i] = k;
            }
        }
        // note: stay matrix doesn't go into maxScore
        // just means we don't end on a stay - oh well...

        // ok, now first check transitions from main matrix
        for (int k=0; k<4; k++)
        {
            if (liks[k] > curlik[i])
            {
                curlik[i] = liks[k];
                // gotta do this to account for implicit zeros
                curstep[i] = likbp[k];
            }
        }

        // now check "transition" back from stay matrix
        if (curstay[i] > curlik[i])
        {
            curlik[i] = curstay[i];
            curstep[i] = U_STAY;
        }

        // and save, if it's the best
        curscore->maxScore <<= MaxInfo(curlik[i],i,refind);
    }
    
    return curscore->maxScore.score;
}

// calculate the score for a given mutation
double Alignment::scoreMutation(const MutInfo& mut, Sequence& mutseq)
{
    // uuhhhh this code is pretty hacky and should not see the
    // light of a yellow sun
	
	// do nothing, no score
	if (this->stripe_width == 0)
		return 0;
    
    // ok, so what it does, is it creates a new column at the back of the
    // scores vector, corresponding to the column at the start of the
    // mutation, and since fillColumn() fills the column _after_ the 
    // backmost column of scores, it'll continue filling the column
    // from there. because we use the column's internal "col" parameter
    // to access the sequence index, instead of its actual index in scores,
    // it more or less works out. but it is a hack, and should be fixed
    
    // save the original size of the vector
    int origsize = scores.size();
    Sequence* origseq = this->sequence;
    // and original score
    //double oldscore = this->getMax();
    double oldscore = columnMax(max(mut.start-3,1));

    this->stripe_width = params.scoring_width;
        
    // and start changing stuff
    this->sequence = &mutseq;
    // add the starting state on to the back
    int startind = max(mut.start-4,0);
    scores.push_back(scores[startind]);
    
    // and fill a few phony columns
    fillColumns(mut.mut.size() + 6);
    
    // and calculate the difference in score. oh god this hurts
    int refind = mut.start+mut.mut.size()+1; // this is the fake forward ind that we need to find
    int fwdind = scores.size()-1;
    while (scores[fwdind]->col > refind && fwdind >= 0)
        fwdind--;
    
    // now scores[fwdind]->col is either equal to refind, or it started off 
    // less than refind, because we were close to the end and didn't 
    // actually fill in enough columns
    // (in which case, we just set refind more aptly)
    if (scores[fwdind]->col >= scores[startind]->col)
        refind = scores[fwdind]->col;
    
    int backind = sequence->states.size()-refind+1;

    double newscore = oldscore-1;
    if (scores[fwdind]->col == refind && fwdind > origsize-1)   
        newscore = columnMax(fwdind,backind);
    else if (params.verbose)
        cerr << "Index mismatch error, mut start" << mut.start << " start " << startind << ", orig " << origsize << ": " << scores[fwdind]->col << " @ " << fwdind << " vs " << refind << ", " << scores.back()->col << " at end " << endl;

    // now roll back to original one
    while (scores.size() > origsize)
        scores.pop_back();
    this->sequence = origseq;

    this->stripe_width = params.realign_width;
    
    // oh gosh, wow....
    return newscore-oldscore;
}

// backtrace the best path to populate ref_align and ref_like
// only needs to be done when matrix is re-populated, not when scoring muts
EventData* Alignment::backtrace()
{
    // note that this vector will have one value for each event state
    // which is the (one-based!) reference sequence index it aligns with
    // (for stays, we just write the same number multiple times)
    vector<int> inds_i;
    vector<int> inds_j;
    vector<double> ref_like;
	
	// do nothing
	if (this->stripe_width == 0)
		return event;
        
    int n0 = event->length;

    int i = scores.back()->maxScore.i;
    int j = scores.back()->maxScore.j;
    int arr = 0; // which array we're on, main or stay

    double score = 1.0;

    while (i>0)
    {
        uint8_t st = *(scores[j]->getStep(i,arr));
        score = *(scores[j]->getPointer(i,arr));

        if (score <= 0.0)
            break;

        switch (st)
        {
            case L_SKIP:
                // skip, don't save
                j--;
                break;
            case UL_MATCH:
                // match, save
                inds_i.push_back(i);
                inds_j.push_back(j);
                ref_like.push_back(score);
                i--;
                j--;
                break;
            case UL_IGNORE:
                // "mismatch", set as insertion
                inds_i.push_back(i);
                inds_j.push_back(-1);
                ref_like.push_back(score);
                i--;
                j--;
                break;                
            case U_INSERT:
                // insert, save but -1 on ref
                inds_i.push_back(i);
                inds_j.push_back(-1);
                ref_like.push_back(score);
                i--;
                break;
            case U_STAY:
                // stay, we're jumping between arrays
                // but only write stuff and change i if it's a stay while
                // in array 1; "stay" in array 0 means we jump to array 1
                // without moving
                if (arr == 1)
                {
                    inds_i.push_back(i);
                    inds_j.push_back(j);
                    ref_like.push_back(score);
                    i--;
                }
                arr = 1-arr;
                break;
            case U_EXTEND:
                // extending stay, save
                inds_i.push_back(i);
                inds_j.push_back(j);
                ref_like.push_back(score);
                i--;
                break;
            case Z_IMPLICIT:
                // we came from an implicit, nonexistant 0
                // (which means skip or match)
                // this is pretty rare, so just quit and leave here
                i = 0;
                break;
            default:
                i = 0;
                break;
        }
    }
    
    // and now refill the internal event with the backtrace results

    // reset arrays
    for (int i=0; i<event->length; i++)
        event->ref_align[i] = 0;
    event->ref_like = event->ref_align;
    // now save the ref values from the backtrace
    for (int i=0; i<inds_i.size(); i++)
    {
        event->ref_align[inds_i[i]-1] = inds_j[i];
        event->ref_like[inds_i[i]-1] = ref_like[i];
    }
    
    // and update the ref_index and other ref stats
    event->updaterefs();
    
    return event;
}
