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

#ifndef _YINUTIL_H_
#define _YINUTIL_H_

#include "vamp-sdk/FFT.h"
#include "MeanFilter.h"
#include "base/Window.h"
#include "maths/MathUtilities.h"

#include <cmath>

#include <iostream>
#include <vector>
#include <exception>

using std::vector;
using std::pair;

class YinUtil
{
public:
    static double sumSquare(const double *in, const size_t startInd, const size_t endInd);
    static void difference(const double *in, double *yinBuffer, const size_t yinBufferSize);
    static void fastDifference(const double *in, double *yinBuffer, const size_t yinBufferSize);
    static void slowDifference(const double *in, double *yinBuffer, const size_t yinBufferSize);
    static void cumulativeDifference(double *yinBuffer, const size_t yinBufferSize);
    static int absoluteThreshold(const double *yinBuffer, const size_t yinBufferSize, const double thresh);
    static vector<double> yinProb(const double *yinBuffer, const size_t prior, const size_t yinBufferSize, size_t minTau = 0, size_t maxTau = 0);
    static double parabolicInterpolation(const double *yinBuffer, const size_t tau,
                                         const size_t yinBufferSize);
    static vector<double> vibrato(const vector<double> pitch, const float inputSampleRate, const size_t stepSize);
    
    struct PeakDetectionParams {
        double _minPos;
        double _maxPos;
        double _threshold;
        double _range;
        bool _interpolate;
        std::string _orderBy;
    };
    
    static vector<pair<double, double> > peakDetection(const vector<double> array,
                                                       PeakDetectionParams params);

private:
    
    static vector<double> movingAverage(vector<double> array, size_t span);
    static void interpolate(const double leftVal, const double middleVal,
                            const double rightVal, int currentBin,
                            double& resultVal, double& resultBin);

    
};

// comparing by magnitude, by default sorts by descending magnitude and in case
// the magnitudes are equal it sorts by ascending position
template<typename Comp1=std::greater<double>,
typename Comp2=std::less_equal<double> >
class ComparePeakMagnitude : public std::binary_function<double, double, bool> {
    Comp1 _cmp1;
    Comp2 _cmp2;
public:
    bool operator () (const pair<double, double>& p1, const pair<double,double>& p2) const {
        if (_cmp1(p1.second, p2.second)) return true;
        if (_cmp1(p2.second, p1.second)) return false;
        return _cmp2(p1.first, p2.first);
    }
};

#endif
