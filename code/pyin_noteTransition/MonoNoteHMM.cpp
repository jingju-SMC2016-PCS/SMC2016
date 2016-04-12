/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    pYIN - A fundamental frequency estimator for monophonic audio
    Centre for Digital Music, Queen Mary, University of London.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "MonoNoteHMM.h"

#include <boost/math/distributions.hpp>

#include <cstdio>
#include <cmath>

using std::vector;
using std::pair;

MonoNoteHMM::MonoNoteHMM() :
    par()
{
    build();
}

const vector<double>
MonoNoteHMM::calculateObsProb(const vector<pair<double, double> > pitchProb)
{
    // pitchProb is a list of pairs (pitches and their probabilities)
    
    size_t nCandidate = pitchProb.size();
    
    // what is the probability of pitched
    double pIsPitched = 0;
    for (size_t iCandidate = 0; iCandidate < nCandidate; ++iCandidate)
    {
        // pIsPitched = pitchProb[iCandidate].second > pIsPitched ? pitchProb[iCandidate].second : pIsPitched;
        pIsPitched += pitchProb[iCandidate].second;
    }

    // pIsPitched = std::pow(pIsPitched, (1-par.priorWeight)) * std::pow(par.priorPitchedProb, par.priorWeight);
    pIsPitched = pIsPitched * (1-par.priorWeight) + par.priorPitchedProb * par.priorWeight;

    vector<double> out = vector<double>(par.n);
    double tempProbSum = 0;
    for (size_t i = 0; i < par.n; ++i)
    {
        if (i % par.nSPP != 2)
        {
            // std::cerr << getMidiPitch(i) << std::endl;
            double tempProb = 0;
            if (nCandidate > 0)
            {
                double minDist = 10000.0;
                double minDistProb = 0;
                size_t minDistCandidate = 0;
                for (size_t iCandidate = 0; iCandidate < nCandidate; ++iCandidate)
                {
                    double currDist = std::abs(getMidiPitch(i)-pitchProb[iCandidate].first);
                    // search minimum distance candidate
                    if (currDist < minDist)
                    {
                        minDist = currDist;
                        minDistProb = pitchProb[iCandidate].second;
                        minDistCandidate = iCandidate;
                    }
                }
                // todo, tempProb * histogram
                tempProb = std::pow(minDistProb, par.yinTrust) * 
                           boost::math::pdf(pitchDistr[i], 
                                            pitchProb[minDistCandidate].first);
            } else {
                tempProb = 1;
            }
            tempProbSum += tempProb;
            out[i] = tempProb;
        }
    }
    
    for (size_t i = 0; i < par.n; ++i)
    {
        if (i % par.nSPP != 2)
        {
            if (tempProbSum > 0) 
            {
                out[i] = out[i] / tempProbSum * pIsPitched;
            }
        } else {
            out[i] = (1-pIsPitched) / (par.nPPS * par.nS);
        }
    }

    return(out);
}

void
MonoNoteHMM::build()
{
    // the states are organised as follows:
    // 0-2. lowest pitch
    //    0. attack state
    //    1. stable state
    //    2. silent state
    // 3-5. second-lowest pitch
    //    3. attack state
    //    ...
    
    // observation distributions
    for (size_t iState = 0; iState < par.n; ++iState)
    {
        pitchDistr.push_back(boost::math::normal(0,1));
        if (iState % par.nSPP == 2)
        {
            // silent state starts tracking
            init.push_back(1.0/(par.nS * par.nPPS));
        } else {
            init.push_back(0.0);            
        }
    }

    for (size_t iPitch = 0; iPitch < (par.nS * par.nPPS); ++iPitch)
    {
        size_t index = iPitch * par.nSPP;
        double mu = par.minPitch + iPitch * 1.0/par.nPPS;
        pitchDistr[index] = boost::math::normal(mu, par.sigmaYinPitchAttack);
        pitchDistr[index+1] = boost::math::normal(mu, par.sigmaYinPitchStable);
        pitchDistr[index+2] = boost::math::normal(mu, 1.0); // dummy
    }
    
    // normal note transition distribution
//    boost::math::normal noteDistanceDistr(0, par.sigma2Note);
    
    static float noteDistanceDistrIndex[115] = {-19.,         -18.66666667, -18.33333333, -18.,         -17.66666667,
                                                -17.33333333, -17.,         -16.66666667, -16.33333333, -16.,         -15.66666667,
                                                -15.33333333, -15.,         -14.66666667, -14.33333333, -14.,         -13.66666667,
                                                -13.33333333, -13.,         -12.66666667, -12.33333333, -12.,         -11.66666667,
                                                -11.33333333, -11.,         -10.66666667, -10.33333333, -10.,          -9.66666667,
                                                -9.33333333,  -9.,          -8.66666667,  -8.33333333,  -8.,          -7.66666667,
                                                -7.33333333,  -7.,          -6.66666667,  -6.33333333,  -6.,          -5.66666667,
                                                -5.33333333,  -5.,          -4.66666667,  -4.33333333,  -4.,          -3.66666667,
                                                -3.33333333,  -3.,          -2.66666667,  -2.33333333,  -2.,          -1.66666667,
                                                -1.33333333,  -1.,          -0.66666667,  -0.33333333,   0.,           0.33333333,
                                                0.66666667,   1.,           1.33333333,   1.66666667,   2.,           2.33333333,
                                                2.66666667,   3.,           3.33333333,   3.66666667,   4.,           4.33333333,
                                                4.66666667,   5.,           5.33333333,   5.66666667,   6.,           6.33333333,
                                                6.66666667,   7.,           7.33333333,   7.66666667,   8.,           8.33333333,
                                                8.66666667,   9.,           9.33333333,   9.66666667,  10.,          10.33333333,
                                                10.66666667,  11.,          11.33333333,  11.66666667, 12.,          12.33333333,
                                                12.66666667,  13.,          13.33333333,  13.66666667,  14.,          14.33333333,
                                                14.66666667,  15.,          15.33333333,  15.66666667,  16.,          16.33333333,
                                                16.66666667,  17.,          17.33333333,  17.66666667,  18.,          18.33333333,
        18.66666667,  19.        };
    
    static float noteDistanceDistr[115] = {  4.80145964e-05,   6.00182455e-05,   7.20218947e-05,   8.40255438e-05,
                                           9.60291929e-05,   1.08032842e-04,   1.20036491e-04,   1.32040140e-04,
                                           1.44043789e-04,   1.56047439e-04,   1.68051088e-04,   1.80054737e-04,
                                           1.92058386e-04,   1.76053520e-04,   1.60048655e-04,   1.44043789e-04,
                                           1.60048654e-04,   1.76053520e-04,   1.92058386e-04,   2.08063251e-04,
                                           2.24068116e-04,   2.40072982e-04,   7.36223812e-04,   1.23237464e-03,
                                           1.72852547e-03,   2.22467630e-03,   2.72082713e-03,   3.21697796e-03,
                                           2.57678334e-03,   1.93658872e-03,   1.29639410e-03,   3.42504121e-03,
                                           5.55368832e-03,   7.68233543e-03,   8.43456411e-03,   9.18679278e-03,
                                           9.93902146e-03,   6.73804837e-03,   3.53707527e-03,   3.36102175e-04,
                                           1.73492742e-02,   3.43624462e-02,   5.13756182e-02,   4.01241978e-02,
                                           2.88727773e-02,   1.76213569e-02,   4.78705526e-02,   7.81197483e-02,
                                           1.08368944e-01,   1.42731390e-01,   1.77093837e-01,   2.11456283e-01,
                                           1.47388806e-01,   8.33213298e-02,   1.92538532e-02,   6.75885468e-02,
                                           1.15923240e-01,   1.64257934e-01,   1.11922024e-01,   5.95861140e-02,
                                           7.25020406e-03,   6.71084010e-02,   1.26966598e-01,   1.86824795e-01,
                                           1.63905828e-01,   1.40986860e-01,   1.18067893e-01,   8.53539478e-02,
                                           5.26400027e-02,   1.99260575e-02,   2.72882956e-02,   3.46505338e-02,
                                           4.20127719e-02,   2.80245195e-02,   1.40362670e-02,   4.80145964e-05,
                                           4.86547910e-03,   9.68294360e-03,   1.45004081e-02,   1.28999216e-02,
                                           1.12994350e-02,   9.69894848e-03,   6.94611162e-03,   4.19327475e-03,
                                           1.44043789e-03,   2.14465197e-03,   2.84886606e-03,   3.55308014e-03,
                                           2.40072982e-03,   1.24837951e-03,   9.60291929e-05,   3.36102175e-04,
                                           5.76175157e-04,   8.16248139e-04,   6.96211648e-04,   5.76175157e-04,
                                           4.56138666e-04,   3.36102175e-04,   2.16065684e-04,   9.60291929e-05,
                                           9.60291929e-05,   9.60291929e-05,   9.60291929e-05,   9.60291929e-05,
                                           9.60291929e-05,   9.60291929e-05,   9.60291929e-05,   9.60291929e-05,
                                           9.60291929e-05,   8.80267602e-05,   8.00243274e-05,   7.20218946e-05,
        6.40194619e-05,   5.60170292e-05,   4.80145964e-05};

    for (size_t iPitch = 0; iPitch < (par.nS * par.nPPS); ++iPitch)
    {
        // loop through all notes and set sparse transition probabilities
        size_t index = iPitch * par.nSPP;

        // transitions from attack state
        from.push_back(index);
        to.push_back(index);
        transProb.push_back(par.pAttackSelftrans);

        from.push_back(index);
        to.push_back(index+1);
        transProb.push_back(1-par.pAttackSelftrans);

        // transitions from stable state
        from.push_back(index+1);
        to.push_back(index+1); // to itself
        transProb.push_back(par.pStableSelftrans);
        
        from.push_back(index+1);
        to.push_back(index+2); // to silent
        transProb.push_back(par.pStable2Silent);

        // the "easy" transitions from silent state
        from.push_back(index+2);
        to.push_back(index+2);
        transProb.push_back(par.pSilentSelftrans);
        
        
        // the more complicated transitions from the silent
        double probSumSilent = 0;

        vector<double> tempTransProbSilent;
        for (size_t jPitch = 0; jPitch < (par.nS * par.nPPS); ++jPitch)
        {
            int fromPitch = iPitch;
            int toPitch = jPitch;
            double semitoneDistance = 
                std::abs(fromPitch - toPitch) * 1.0 / par.nPPS;
            
            // if (std::fmod(semitoneDistance, 1) == 0 && semitoneDistance > par.minSemitoneDistance)
            if (semitoneDistance == 0 || 
                (semitoneDistance > par.minSemitoneDistance 
                 && semitoneDistance < par.maxJump))
            {
                size_t toIndex = jPitch * par.nSPP; // note attack index

//                double tempWeightSilent = boost::math::pdf(noteDistanceDistr,
//                                                           semitoneDistance);
                
                size_t semitoneDistanceIndex = 0;
                for (size_t iNDDI = 1; iNDDI < sizeof(noteDistanceDistrIndex)-1 ; ++iNDDI) {
                    float diffPrevious    = std::abs(noteDistanceDistrIndex[iNDDI-1]  - semitoneDistance);
                    float diffCurrent     = std::abs(noteDistanceDistrIndex[iNDDI]    - semitoneDistance);
                    float diffNext        = std::abs(noteDistanceDistrIndex[iNDDI+1]  - semitoneDistance);
                    
                    if (diffCurrent < diffPrevious && diffCurrent < diffNext) {
                        semitoneDistanceIndex = iNDDI;
                        break;
                    }
                }
                
                double tempWeightSilent = noteDistanceDistr[semitoneDistanceIndex];
                
                // debug
                //std::cerr << semitoneDistance << " -- "<< semitoneDistanceIndex << " -- "<< tempWeightSilent << std::endl;
                
                probSumSilent += tempWeightSilent;

                tempTransProbSilent.push_back(tempWeightSilent);

                from.push_back(index+2);
                to.push_back(toIndex);
            }
        }
        for (size_t i = 0; i < tempTransProbSilent.size(); ++i)
        {
            transProb.push_back((1-par.pSilentSelftrans) * tempTransProbSilent[i]/probSumSilent);
        }
    }
}

double
MonoNoteHMM::getMidiPitch(size_t index)
{
    return pitchDistr[index].mean();
}

double
MonoNoteHMM::getFrequency(size_t index)
{
    return 440 * pow(2, (pitchDistr[index].mean()-69)/12);
}
