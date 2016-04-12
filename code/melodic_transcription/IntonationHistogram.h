/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
 * Copyright (C) 2016  Music Technology Group - Universitat Pompeu Fabra
 *
 * This file is part of pYIN Beijing Opera version
 *
 * pYIN Beijing Opera version is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the Free
 * Software Foundation (FSF), either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the Affero GNU General Public License
 * version 3 along with this program.  If not, see http://www.gnu.org/licenses/
 *
 * If you have any problem about this version code, please contact: Rong Gong
 * rong.gong@upf.edu
 *
 * If you have any problem about this algorithm, I suggest you to contact: Matthias Mauch
 * m.mauch@qmul.ac.uk who is the original C++ version author of this algorithm
 *
 * If you want to refer this code, please consider this article:
 *
 * M. Mauch and S. Dixon,
 * “pYIN: A Fundamental Frequency Estimator Using Probabilistic Threshold Distributions”,
 * in Proceedings of the IEEE International Conference on Acoustics,
 * Speech, and Signal Processing (ICASSP 2014), 2014.
 *
 * M. Mauch, C. Cannam, R. Bittner, G. Fazekas, J. Salamon, J. Dai, J. Bello and S. Dixon,
 * “Computer-aided Melody Note Transcription Using the Tony Software: Accuracy and Efficiency”,
 * in Proceedings of the First International Conference on Technologies for
 * Music Notation and Representation, 2015.
 
 */

#ifndef _INTONATIONHISTOGRAM_H_
#define _INTONATIONHISTOGRAM_H_


#include <vector>
#include <cstdio>
#include <iostream>

using std::vector;
using std::pair;

class IntonationHistogram
{
public:
    IntonationHistogram();
    const std::vector<std::pair<double,double> > histogram(const vector<double> pitch);//, const vector<double> level);
    const std::vector<std::pair<double,double> > histogramFusion(vector<pair<double,double> > mpHisto, vector<pair<double,double> > npHisto, float fusionWeightStd);
    // double getMidiPitch(size_t index);
    // double getFrequency(size_t index);
    double m_minFreq; // 82.40689f/2
    size_t m_nBPS;
    size_t m_nPitch;
    vector<double> m_freqs;
};


#endif

