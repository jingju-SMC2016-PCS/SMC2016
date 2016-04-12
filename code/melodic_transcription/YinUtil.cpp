/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    pYIN - A fundamental frequency estimator for monophonic audio
    Centre for Digital Music, Queen Mary, University of London.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
 
    last update 03/03/2016: Rong Gong, MTG-UPF, rong.gong@upf.edu
    add moving average, peakDetection and vibrato section detection
*/

#include "YinUtil.h"

#include <vector>

#include <cstdio>
#include <cmath>
#include <algorithm>

#include <boost/math/distributions.hpp>

void 
YinUtil::slowDifference(const double *in, double *yinBuffer, const size_t yinBufferSize) 
{
    yinBuffer[0] = 0;
    double delta ;
    int startPoint = 0;
    int endPoint = 0;
    for (int i = 1; i < yinBufferSize; ++i) {
        yinBuffer[i] = 0;
        startPoint = yinBufferSize/2 - i/2;
        endPoint = startPoint + yinBufferSize;
        for (int j = startPoint; j < endPoint; ++j) {
            delta = in[i+j] - in[j];
            yinBuffer[i] += delta * delta;
        }
    }    
}

void 
YinUtil::fastDifference(const double *in, double *yinBuffer, const size_t yinBufferSize) 
{
    
    // DECLARE AND INITIALISE
    // initialisation of most of the arrays here was done in a separate function,
    // with all the arrays as members of the class... moved them back here.
    
    size_t frameSize = 2 * yinBufferSize;
    
    double *audioTransformedReal = new double[frameSize];
    double *audioTransformedImag = new double[frameSize];
    double *nullImag = new double[frameSize];
    double *kernel = new double[frameSize];
    double *kernelTransformedReal = new double[frameSize];
    double *kernelTransformedImag = new double[frameSize];
    double *yinStyleACFReal = new double[frameSize];
    double *yinStyleACFImag = new double[frameSize];
    double *powerTerms = new double[yinBufferSize];
    
    for (size_t j = 0; j < yinBufferSize; ++j)
    {
        yinBuffer[j] = 0.; // set to zero
        powerTerms[j] = 0.; // set to zero
    }
    
    for (size_t j = 0; j < frameSize; ++j)
    {
        nullImag[j] = 0.;
        audioTransformedReal[j] = 0.;
        audioTransformedImag[j] = 0.;
        kernel[j] = 0.;
        kernelTransformedReal[j] = 0.;
        kernelTransformedImag[j] = 0.;
        yinStyleACFReal[j] = 0.;
        yinStyleACFImag[j] = 0.;
    }
    
    // POWER TERM CALCULATION
    // ... for the power terms in equation (7) in the Yin paper
    powerTerms[0] = 0.0;
    for (size_t j = 0; j < yinBufferSize; ++j) {
        powerTerms[0] += in[j] * in[j];
    }

    // now iteratively calculate all others (saves a few multiplications)
    for (size_t tau = 1; tau < yinBufferSize; ++tau) {
        powerTerms[tau] = powerTerms[tau-1] - in[tau-1] * in[tau-1] + in[tau+yinBufferSize] * in[tau+yinBufferSize];  
    }
    
    // YIN-STYLE AUTOCORRELATION via FFT
    // 1. data
    Vamp::FFT::forward(frameSize, in, nullImag, audioTransformedReal, audioTransformedImag);
    
    // 2. half of the data, disguised as a convolution kernel
    for (size_t j = 0; j < yinBufferSize; ++j) {
        kernel[j] = in[yinBufferSize-1-j];
    }
    Vamp::FFT::forward(frameSize, kernel, nullImag, kernelTransformedReal, kernelTransformedImag);

    // 3. convolution via complex multiplication -- written into
    for (size_t j = 0; j < frameSize; ++j) {
        yinStyleACFReal[j] = audioTransformedReal[j]*kernelTransformedReal[j] - audioTransformedImag[j]*kernelTransformedImag[j]; // real
        yinStyleACFImag[j] = audioTransformedReal[j]*kernelTransformedImag[j] + audioTransformedImag[j]*kernelTransformedReal[j]; // imaginary
    }
    Vamp::FFT::inverse(frameSize, yinStyleACFReal, yinStyleACFImag, audioTransformedReal, audioTransformedImag);
    
    // CALCULATION OF difference function
    // ... according to (7) in the Yin paper.
    for (size_t j = 0; j < yinBufferSize; ++j) {
        // taking only the real part
        yinBuffer[j] = powerTerms[0] + powerTerms[j] - 2 * audioTransformedReal[j+yinBufferSize-1];
    }
    delete [] audioTransformedReal;
    delete [] audioTransformedImag;
    delete [] nullImag;
    delete [] kernel;
    delete [] kernelTransformedReal;
    delete [] kernelTransformedImag;
    delete [] yinStyleACFReal;
    delete [] yinStyleACFImag;
    delete [] powerTerms;
}

void 
YinUtil::cumulativeDifference(double *yinBuffer, const size_t yinBufferSize)
{    
    size_t tau;
    
    yinBuffer[0] = 1;
    
    double runningSum = 0;
    
    for (tau = 1; tau < yinBufferSize; ++tau) {
        runningSum += yinBuffer[tau];
        if (runningSum == 0)
        {
            yinBuffer[tau] = 1;
        } else {
            yinBuffer[tau] *= tau / runningSum;
        }
    }    
}

int 
YinUtil::absoluteThreshold(const double *yinBuffer, const size_t yinBufferSize, const double thresh)
{
    size_t tau;
    size_t minTau = 0;
    double minVal = 1000.;
    
    // using Joren Six's "loop construct" from TarsosDSP
    tau = 2;
    while (tau < yinBufferSize)
    {
        if (yinBuffer[tau] < thresh)
        {
            while (tau+1 < yinBufferSize && yinBuffer[tau+1] < yinBuffer[tau])
            {
                ++tau;
            }
            return tau;
        } else {
            if (yinBuffer[tau] < minVal)
            {
                minVal = yinBuffer[tau];
                minTau = tau;
            }
        }
        ++tau;
    }
    if (minTau > 0)
    {
        return -minTau;
    }
    return 0;
}

static float uniformDist[100] = {0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000,0.0100000};
static float betaDist1[100] = {0.028911,0.048656,0.061306,0.068539,0.071703,0.071877,0.069915,0.066489,0.062117,0.057199,0.052034,0.046844,0.041786,0.036971,0.032470,0.028323,0.024549,0.021153,0.018124,0.015446,0.013096,0.011048,0.009275,0.007750,0.006445,0.005336,0.004397,0.003606,0.002945,0.002394,0.001937,0.001560,0.001250,0.000998,0.000792,0.000626,0.000492,0.000385,0.000300,0.000232,0.000179,0.000137,0.000104,0.000079,0.000060,0.000045,0.000033,0.000024,0.000018,0.000013,0.000009,0.000007,0.000005,0.000003,0.000002,0.000002,0.000001,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000};
static float betaDist2[100] = {0.012614,0.022715,0.030646,0.036712,0.041184,0.044301,0.046277,0.047298,0.047528,0.047110,0.046171,0.044817,0.043144,0.041231,0.039147,0.036950,0.034690,0.032406,0.030133,0.027898,0.025722,0.023624,0.021614,0.019704,0.017900,0.016205,0.014621,0.013148,0.011785,0.010530,0.009377,0.008324,0.007366,0.006497,0.005712,0.005005,0.004372,0.003806,0.003302,0.002855,0.002460,0.002112,0.001806,0.001539,0.001307,0.001105,0.000931,0.000781,0.000652,0.000542,0.000449,0.000370,0.000303,0.000247,0.000201,0.000162,0.000130,0.000104,0.000082,0.000065,0.000051,0.000039,0.000030,0.000023,0.000018,0.000013,0.000010,0.000007,0.000005,0.000004,0.000003,0.000002,0.000001,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000};
static float betaDist3[100] = {0.006715,0.012509,0.017463,0.021655,0.025155,0.028031,0.030344,0.032151,0.033506,0.034458,0.035052,0.035331,0.035332,0.035092,0.034643,0.034015,0.033234,0.032327,0.031314,0.030217,0.029054,0.027841,0.026592,0.025322,0.024042,0.022761,0.021489,0.020234,0.019002,0.017799,0.016630,0.015499,0.014409,0.013362,0.012361,0.011407,0.010500,0.009641,0.008830,0.008067,0.007351,0.006681,0.006056,0.005475,0.004936,0.004437,0.003978,0.003555,0.003168,0.002814,0.002492,0.002199,0.001934,0.001695,0.001481,0.001288,0.001116,0.000963,0.000828,0.000708,0.000603,0.000511,0.000431,0.000361,0.000301,0.000250,0.000206,0.000168,0.000137,0.000110,0.000088,0.000070,0.000055,0.000043,0.000033,0.000025,0.000019,0.000014,0.000010,0.000007,0.000005,0.000004,0.000002,0.000002,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000};
static float betaDist4[100] = {0.003996,0.007596,0.010824,0.013703,0.016255,0.018501,0.020460,0.022153,0.023597,0.024809,0.025807,0.026607,0.027223,0.027671,0.027963,0.028114,0.028135,0.028038,0.027834,0.027535,0.027149,0.026687,0.026157,0.025567,0.024926,0.024240,0.023517,0.022763,0.021983,0.021184,0.020371,0.019548,0.018719,0.017890,0.017062,0.016241,0.015428,0.014627,0.013839,0.013068,0.012315,0.011582,0.010870,0.010181,0.009515,0.008874,0.008258,0.007668,0.007103,0.006565,0.006053,0.005567,0.005107,0.004673,0.004264,0.003880,0.003521,0.003185,0.002872,0.002581,0.002312,0.002064,0.001835,0.001626,0.001434,0.001260,0.001102,0.000959,0.000830,0.000715,0.000612,0.000521,0.000440,0.000369,0.000308,0.000254,0.000208,0.000169,0.000136,0.000108,0.000084,0.000065,0.000050,0.000037,0.000027,0.000019,0.000014,0.000009,0.000006,0.000004,0.000002,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000};
static float single10[100] = {0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,1.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000};
static float single15[100] = {0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,1.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000};
static float single20[100] = {0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,1.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000};

std::vector<double>
YinUtil::yinProb(const double *yinBuffer, const size_t prior, const size_t yinBufferSize, const size_t minTau0, const size_t maxTau0) 
{
    size_t minTau = 2;
    size_t maxTau = yinBufferSize;

    // adapt period range, if necessary
    if (minTau0 > 0 && minTau0 < maxTau0) minTau = minTau0;
    if (maxTau0 > 0 && maxTau0 < yinBufferSize && maxTau0 > minTau) maxTau = maxTau0;

    double minWeight = 0.01;
    size_t tau;
    std::vector<float> thresholds;
    std::vector<float> distribution;
    std::vector<double> peakProb = std::vector<double>(yinBufferSize);
    
    size_t nThreshold = 100;
    int nThresholdInt = nThreshold;
    
    for (int i = 0; i < nThresholdInt; ++i)
    {
        switch (prior) {
            case 0:
                distribution.push_back(uniformDist[i]);
                break;
            case 1:
                distribution.push_back(betaDist1[i]);
                break;
            case 2:
                distribution.push_back(betaDist2[i]);
                break;
            case 3:
                distribution.push_back(betaDist3[i]);
                break;
            case 4:
                distribution.push_back(betaDist4[i]);
                break;
            case 5:
                distribution.push_back(single10[i]);
                break;
            case 6:
                distribution.push_back(single15[i]);
                break;
            case 7:
                distribution.push_back(single20[i]);
                break;
            default:
                distribution.push_back(uniformDist[i]);
        }
        thresholds.push_back(0.01 + i*0.01);
    }
    
    
    int currThreshInd = nThreshold-1;
    tau = minTau;
    
    // double factor = 1.0 / (0.25 * (nThresholdInt+1) * (nThresholdInt + 1)); // factor to scale down triangular weight
    size_t minInd = 0;
    float minVal = 42.f;
    // while (currThreshInd != -1 && tau < maxTau)
    // {
    //     if (yinBuffer[tau] < thresholds[currThreshInd])
    //     {
    //         while (tau + 1 < maxTau && yinBuffer[tau+1] < yinBuffer[tau])
    //         {
    //             tau++;
    //         }
    //         // tau is now local minimum
    //         // std::cerr << tau << " " << currThreshInd << " "<< thresholds[currThreshInd] << " " << distribution[currThreshInd] << std::endl;
    //         if (yinBuffer[tau] < minVal && tau > 2){
    //             minVal = yinBuffer[tau];
    //             minInd = tau;
    //         }
    //         peakProb[tau] += distribution[currThreshInd];
    //         currThreshInd--;
    //     } else {
    //         tau++;
    //     }
    // }
    // double nonPeakProb = 1;
    // for (size_t i = minTau; i < maxTau; ++i)
    // {
    //     nonPeakProb -= peakProb[i];
    // }
    // 
    // std::cerr << tau << " " << currThreshInd << " "<< thresholds[currThreshInd] << " " << distribution[currThreshInd] << std::endl;
    float sumProb = 0;
    while (tau+1 < maxTau)
    {
        if (yinBuffer[tau] < thresholds[thresholds.size()-1] && yinBuffer[tau+1] < yinBuffer[tau])
        {
            while (tau + 1 < maxTau && yinBuffer[tau+1] < yinBuffer[tau])
            {
                tau++;
            }
            // tau is now local minimum
            // std::cerr << tau << " " << currThreshInd << " "<< thresholds[currThreshInd] << " " << distribution[currThreshInd] << std::endl;
            if (yinBuffer[tau] < minVal && tau > 2){
                minVal = yinBuffer[tau];
                minInd = tau;
            }
            currThreshInd = nThresholdInt-1;
            while (thresholds[currThreshInd] > yinBuffer[tau] && currThreshInd > -1) {
                // std::cerr << distribution[currThreshInd] << std::endl;
                peakProb[tau] += distribution[currThreshInd];
                currThreshInd--;
            }
            // peakProb[tau] = 1 - yinBuffer[tau];
            sumProb += peakProb[tau];
            tau++;
        } else {
            tau++;
        }
    }
    
    if (peakProb[minInd] > 1) {
        std::cerr << "WARNING: yin has prob > 1 ??? I'm returning all zeros instead." << std::endl;
        return(std::vector<double>(yinBufferSize));
    }
    
    double nonPeakProb = 1;
    if (sumProb > 0) {
        for (size_t i = minTau; i < maxTau; ++i)
        {
            peakProb[i] = peakProb[i] / sumProb * peakProb[minInd];
            nonPeakProb -= peakProb[i];
        }
    }
    if (minInd > 0)
    {
        // std::cerr << "min set " << minVal << " " << minInd << " " << nonPeakProb << std::endl; 
        peakProb[minInd] += nonPeakProb * minWeight;
    }
    
    return peakProb;
}

double
YinUtil::parabolicInterpolation(const double *yinBuffer, const size_t tau, const size_t yinBufferSize) 
{
    // this is taken almost literally from Joren Six's Java implementation
    if (tau == yinBufferSize) // not valid anyway.
    {
        return static_cast<double>(tau);
    }
    
    double betterTau = 0.0;
    if (tau > 0 && tau < yinBufferSize-1) {
        float s0, s1, s2;
        s0 = yinBuffer[tau-1];
        s1 = yinBuffer[tau];
        s2 = yinBuffer[tau+1];
        
        double adjustment = (s2 - s0) / (2 * (2 * s1 - s2 - s0));
        
        if (abs(adjustment)>1) adjustment = 0;
        
        betterTau = tau + adjustment;
    } else {
        // std::cerr << "WARNING: can't do interpolation at the edge (tau = " << tau << "), will return un-interpolated value.\n";
        betterTau = tau;
    }
    return betterTau;
}

double 
YinUtil::sumSquare(const double *in, const size_t start, const size_t end)
{
    double out = 0;
    for (size_t i = start; i < end; ++i)
    {
        out += in[i] * in[i];
    }
    return out;
}

vector<double>
YinUtil::movingAverage(vector<double> array, size_t span)
{
    /* matlab smooth function */
    assert(span%2 != 0); // span should be odd
    
    size_t halfSpan = size_t(span/2);
    size_t tempHalfSpan;
    double averge;
    vector<double> tempArray;
    vector<double> outArray;
    
    for (size_t iArray = 0; iArray < array.size(); ++iArray) {
        tempArray.clear();
        if (int(iArray)-int(halfSpan) < 0) {
            tempHalfSpan = iArray;
        } else if (iArray+halfSpan >= array.size()) {
            tempHalfSpan = array.size()-iArray-1;
        } else {
            tempHalfSpan = halfSpan;
        }
        
        for (size_t iiArray = iArray - tempHalfSpan; iiArray < iArray + tempHalfSpan + 1; ++iiArray) {
            tempArray.push_back(array[iiArray]);
//            std::cerr << iiArray << " temp Array " << array[iiArray] << std::endl;

        }
        
//        std::cin.get();
        
        averge = std::accumulate(tempArray.begin(), tempArray.end(), 0.0)/tempArray.size();
        outArray.push_back(averge);
    }
    return(outArray);
}

void
YinUtil::interpolate(const double leftVal, const double middleVal, const double rightVal, int currentBin, double& resultVal, double& resultBin)
{
    double delta_x = 0.5 * ((leftVal - rightVal) / (leftVal - 2*middleVal + rightVal));
    resultBin = currentBin + delta_x;
    resultVal = middleVal - 0.25 * (leftVal - rightVal) * delta_x;
}

vector<pair<double, double> >
YinUtil::peakDetection(const vector<double> array, PeakDetectionParams params)
{
    double _minPos      = params._minPos;
    double _maxPos      = params._maxPos;
    double _threshold   = params._threshold;
    double _range       = params._range;
    bool _interpolate   = params._interpolate;
    std::string _orderBy = params._orderBy;
    
    const int size = (int)array.size();
    
    assert(size >= 2); // array size should > 2
    
    const double scale = _range / (double)(size - 1);
    
    vector<pair<double, double> > peaks;
    pair<double, double> peak;
    
    // we want to round up to the next integer instead of simple truncation,
    // otherwise the peak frequency at i can be lower than _minPos
    int i = std::max(0, (int) ceil(_minPos / scale));
    
    // first check the boundaries:
    if (i+1 < size && array[i] > array[i+1]) {
        if (array[i] > _threshold) {
            peak.first  = i*scale;
            peak.second = array[i];
            peaks.push_back(peak);
        }
    }

    while(true) {
        // going down
        while (i+1 < size-1 && array[i] >= array[i+1]) {
            i++;
        }
        
        // now we're climbing
        while (i+1 < size-1 && array[i] < array[i+1]) {
            i++;
        }
        
        // not anymore, go through the plateau
        int j = i;
        while (j+1 < size-1 && (array[j] == array[j+1])) {
            j++;
        }
        
        // end of plateau, do we go up or down?
        if (j+1 < size-1 && array[j+1] < array[j] && array[j] > _threshold) { // going down again
            double resultBin = 0.0;
            double resultVal = 0.0;
            
            if (j != i) { // plateau peak between i and j
                if (_interpolate) {
                    resultBin = (i + j) * 0.5;
                }
                else {
                    resultBin = i;
                }
                resultVal = array[i];
            }
            else { // interpolate peak at i-1, i and i+1
                if (_interpolate) {
                    interpolate(array[j-1], array[j], array[j+1], j, resultVal, resultBin);
                }
                else {
                    resultBin = j;
                    resultVal = array[j];
                }
            }
            
            double resultPos = resultBin * scale;
            
            if (resultPos > _maxPos)
                break;
            
            peak.first  = resultPos;
            peak.second = resultVal;
            peaks.push_back(peak);
        }
        
        // nothing found, start loop again
        i = j;
        
        if (i+1 >= size-1) { // check the one just before the last position
            if (i == size-2 && array[i-1] < array[i] &&
                array[i+1] < array[i] &&
                array[i] > _threshold) {
                double resultBin = 0.0;
                double resultVal = 0.0;
                if (_interpolate) {
                    interpolate(array[i-1], array[i], array[i+1], j, resultVal, resultBin);
                }
                else {
                    resultBin = i;
                    resultVal = array[i];
                }
                
                peak.first  = resultBin*scale;
                peak.second = resultVal;
                peaks.push_back(peak);
            }
            break;
        }
    }
    
    // check upper boundary here, so peaks are already sorted by position
    float pos = _maxPos/scale;
    if (size-2 <pos && pos <= size-1 && array[size-1] > array[size-2]) {
        if (array[size-1] > _threshold) {
            peak.first = (size-1)*scale;
            peak.second = array[size-1];
            peaks.push_back(peak);
        }
    }
    
    if (_orderBy == "amplitude") {
        // sort peaks by magnitude, in case of equality,
        // return the one having smaller position
        std::sort(peaks.begin(), peaks.end(),
                  ComparePeakMagnitude<std::greater<double>, std::less<double> >());
    }
    else if (_orderBy == "position") {
        // they're already sorted by position
    }
    else {
        throw;
    }
    
    return(peaks);
}


std::vector<double>
YinUtil::vibrato(const vector<double> pitch, const float inputSampleRate, const size_t stepSize)
{
//    vector<double> out;
    
    float sampleRate = inputSampleRate/stepSize;
    
    size_t frameSize = size_t(std::roundf(0.5*sampleRate));
    
    if (pitch.size() < frameSize) { //  dynamic frameSize
        frameSize = pitch.size();
    }

    int fftSize = MathUtilities::nextPowerOfTwo(4*frameSize);
    float dBLobe = 15;
    float dBSecondLobe = 20;
    float minFreq = 2.0;
    float maxFreq = 8.0;
    
    // Initialise FFT object
    Window<double> *window;
    double *windowbuf;
    
    double *in          = new double[fftSize];
    double *nullImag    = new double[fftSize];
    double *m_FFTRe     = new double[fftSize];
    double *m_FFTIm     = new double[fftSize];
    vector<double> m_FFTMag;
    
    window = 0;
    windowbuf = 0;
    
    double minExt = 30.0;
//    maxExt = 250.0;
//    fRef = 55.0;
//    winAnalysis = 'hann';
    
    // INIT
    vector<double> f0;
    for (size_t iPitch = 0; iPitch < pitch.size(); ++iPitch) {
        if (pitch[iPitch] < 0) {
            f0.push_back(0);
        } else {
            f0.push_back(pitch[iPitch]);
        }
    }
    
    // GET CONTOUR SEGMENTS
    vector<size_t> startC;
    vector<size_t> endC;
    
    if (f0[0]>0) {
        startC.push_back(0);
    }
    
    for (size_t iF0 = 0; iF0 < f0.size()-1; ++iF0) {
        if (std::abs(f0[iF0+1]) > 0 && f0[iF0] == 0) {
            startC.push_back(iF0+1);
        }
        if (f0[iF0+1]==0 && std::abs(f0[iF0])>0) {
            endC.push_back(iF0);
        }
    }
    
    if (endC.size() < startC.size()) {
        endC.push_back(f0.size());
    }
    
    // vibrato annotations
    vector<double> vibSec = vector<double>(f0.size());
    
//    vibFreq = [0]*len(f0)
//    vibFreqMin = [0]*len(f0)
//    vibFreqMax = [0]*len(f0)
//    vibExt = [0]*len(f0)
    
    // ANALYSE EACH SEGMENT
    for (size_t iC = 0; iC < startC.size(); ++iC) {
        
        // get segment in cents
        vector<double> contour(f0.begin()+startC[iC], f0.begin()+endC[iC]);
        
        // frame-wise FFT
        size_t iFrame = 0;

        while (iFrame <= contour.size() - frameSize) {
            
            m_FFTMag.resize(0);
            for (size_t iiFrame = 0; iiFrame < fftSize; ++iiFrame) {
                in[iiFrame]         = 0.0;
                nullImag[iiFrame]   = 0.0;
                m_FFTRe[iiFrame]    = 0.0;
                m_FFTIm[iiFrame]    = 0.0;
            }
            
            vector<double> frame(contour.begin()+iFrame, contour.begin()+iFrame+frameSize);
            
            // extend constraint
            double maxFrame = std::pow(2,(*std::max_element(frame.begin(), frame.end()) - 69) / 12) * 440;
            double minFrame = std::pow(2,(*std::min_element(frame.begin(), frame.end()) - 69) / 12) * 440;
            
            double extent = maxFrame - minFrame;
            
//            std::cerr << " extent " << extent << " " << frame.size() << std::endl;
            
            if (extent < minExt) {
                iFrame += size_t(std::roundf(frameSize/2.0));
                continue;
            }
            
            // remove dc
            double meanFrame = MathUtilities::mean(frame);
            
            for (size_t iiFrame = 0; iiFrame < frame.size(); ++iiFrame) {
                frame[iiFrame] -= meanFrame;
            }
            
            // moving average
            frame = movingAverage(frame, 5);
            
//            for (size_t iiFrame = 0; iiFrame < frame.size(); ++iiFrame) {
//                std::cerr << " moving average " << iiFrame << " " << frame[iiFrame] << " " << startC[iC] << " " << iFrame << " " << contour.size()-frameSize << " " << iFrame <<  std::endl;
//            }
//            
//            std::cin.get();
            
            // windowing
            if (!window) {
                window = new Window<double>(HanningWindow, frameSize);
                windowbuf = new double[frameSize];
            }
            
            for (size_t iiFrame = 0; iiFrame < frame.size(); ++iiFrame) {
                windowbuf[iiFrame] = frame[iiFrame];
            }
            
            window->cut(windowbuf);
            
            // remove again dc
            double meanWindowbuf = MathUtilities::mean(windowbuf, frameSize);
            
            for (size_t iiFrame = 0; iiFrame < frameSize; ++iiFrame) {
                windowbuf[iiFrame] -= meanWindowbuf;
            }

            // zero phrase padding
            for (size_t iiFrame = size_t(frame.size()/2); iiFrame < frame.size(); ++iiFrame) {
                in[iiFrame - size_t(frame.size()/2)] = windowbuf[iiFrame];
            }
            
            for (size_t iiFrame = 0; iiFrame < size_t(frame.size()/2); ++iiFrame) {
                in[fftSize - size_t(frame.size()/2) + iiFrame] = windowbuf[iiFrame];
            }
            
//            for (size_t iiFrame = 0; iiFrame < fftSize; ++iiFrame) {
//                std::cerr << iiFrame << " in " << in[iiFrame] << std::endl;
//            }
            
            // FFT magnitude
            Vamp::FFT::forward(fftSize, in, nullImag, m_FFTRe, m_FFTIm);
            
            for (size_t iiFrame = 0; iiFrame < fftSize/2+1; ++iiFrame) {
                m_FFTMag.push_back(std::sqrt(m_FFTRe[iiFrame] * m_FFTRe[iiFrame]
                                              + m_FFTIm[iiFrame] * m_FFTIm[iiFrame]));
//                std::cerr << iiFrame << " m_FFTMag " << m_FFTMag[iiFrame] << std::endl;
            }
            
//            std::cin.get();
            
            // peak detection parameters
            PeakDetectionParams params;
            params._minPos = 0.;
            params._maxPos = 1.;
            params._threshold = -1e6;
            params._range = 1.0;
            params._interpolate = true;
            params._orderBy = "amplitude";
            
            vector<pair<double, double> > peaks = peakDetection(m_FFTMag, params);
            
//            size_t nWantedPeaks = peaks.size();
            size_t nWantedPeaks = size_t(std::fmin(peaks.size(), 3));
            vector<double> freqs(nWantedPeaks);
            vector<double> amp(nWantedPeaks);
            
            for (size_t iPeak = 0; iPeak < nWantedPeaks; ++iPeak) {
                freqs[iPeak]    = peaks[iPeak].first*(fftSize/2+1)*sampleRate/fftSize;
                amp[iPeak]      = peaks[iPeak].second;
//                std::cerr << iPeak << " peak " << freqs[iPeak] << " " << amp[iPeak] << std::endl;
            }
            
//            std::cin.get();
            
            //print freqs, amp
            if (freqs.size()<=0) {
                iFrame += size_t(std::roundf(frameSize/2.0));
                continue;
            }
            if (freqs[0]<minFreq || freqs[0]>maxFreq) { // strongest peak is not in considered range
                iFrame += size_t(std::roundf(frameSize/2.0));
                continue;
            }
            if (freqs.size()>1) { // there is a second peak
                if (freqs[1]>minFreq && freqs[1]<maxFreq) {
                    iFrame += size_t(std::roundf(frameSize/2.0));
                    continue;
                }
                if (20*std::log10(amp[0]/amp[1])<dBLobe) {
                    iFrame += size_t(std::roundf(frameSize/2.0));
                    continue;
                }
            }
            if (freqs.size()>2) { //there is a third peak
                if (freqs[2]>minFreq && freqs[2]<maxFreq) { // it is also in the vibrato range
                    iFrame += size_t(std::roundf(frameSize/2.0));
                    continue;
                }
                if (20*std::log10(amp[0]/amp[2])<dBSecondLobe) {
                    iFrame += size_t(std::roundf(frameSize/2.0));
                    continue;
                }
            }

            for (size_t iiFrame = startC[iC]+iFrame; iiFrame < startC[iC]+iFrame+frameSize; ++iiFrame) {
                vibSec[iiFrame] = f0[iiFrame];
            }
            
            iFrame += size_t(std::roundf(frameSize/2.0));
            
        }
    }
    
    delete[] windowbuf;
    delete window;
    
    delete[] in;
    delete[] nullImag;
    delete[] m_FFTRe;
    delete[] m_FFTIm;
    
    return(vibSec);
}

