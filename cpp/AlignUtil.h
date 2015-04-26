/*==========================================================
 *
 * AlignUtil.h
 *
 * Helper function definitions for all the various functions.
 *
 *========================================================*/

#ifndef _ALIGNUTIL_H_
#define _ALIGNUTIL_H_

#include <cmath>
#include <stdint.h>
#include <string>

using namespace std;


const int N_STATES = 1024;
const double inf=1e300;
#ifndef M_PI
const double M_PI=3.1415926;
#endif
const double log2pi = log(2*M_PI);

// gaussian dist functions
inline double normpdf(double x, double mu, double sigma)
{
    double d = (x-mu)/sigma;
    d = exp(-0.5*d*d)/sigma/sqrt(2*M_PI);
    return d;
}

inline double lognormpdf(double x, double mu, double sigma, double logsigma)
{
    double d = (x-mu)/sigma;
    return -0.5*(d*d + log2pi) - logsigma;
}

// inverse gaussian dist functions
inline double igpdf(double x, double mu, double lambda)
{
    double d = (x-mu)/mu;
    d = exp(-0.5*d*d*lambda/x);
    return d*sqrt(lambda/(2*M_PI*x*x*x));
}

inline double logigpdf(double x, double mu, double lambda, double logx, double loglambda)
{
    double d = (x-mu)/mu;
    //return 0.5*(log(lambda/(x*x*x)) - log2pi - d*d*lambda/x);
    return 0.5*(loglambda - 3*logx - log2pi - d*d*lambda/x);
}


// Contains global parameters for alignment algorithm
struct AlignParams
{
    double lik_offset;
    int scoring_width;
    int realign_width;
    int verbose;
    
    AlignParams() : lik_offset(4.5), scoring_width(150), realign_width(300), verbose(0)
    {}
};

// little struct to hold mutation information
struct MutInfo
{
    int start;
    string orig;
    string mut;

    MutInfo() : start(0)
    {}
};

// same as above, but with a score
struct MutScore : MutInfo
{
    double score;

    // initialize with a slightly negative score to get rid of null
    // mutations, which are typically ~+-1e-12    
    MutScore(const MutInfo& mut) : MutInfo(mut), score(-1e-6)
    {}

    MutScore() : MutInfo(), score(-1e-6)
{}
};

#endif