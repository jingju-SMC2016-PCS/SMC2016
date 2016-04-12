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

#include "IntonationHistogram.h"

#include <boost/math/distributions.hpp>

#include <cstdio>
#include <cmath>

using std::vector;
using std::pair;

IntonationHistogram::IntonationHistogram() :
m_minFreq(61.735),
m_nBPS(5),
m_nPitch(0),
m_freqs(0)
{
    m_nPitch = 69 * m_nBPS;
    m_freqs = vector<double>(m_nPitch);
    for (size_t iPitch = 0; iPitch < m_nPitch; ++iPitch)
    {
        m_freqs[iPitch] = m_minFreq * std::pow(2, iPitch * 1.0 / (12 * m_nBPS));
        m_freqs[iPitch] = 12 * std::log(m_freqs[iPitch]/440)/std::log(2.) + 69;
        
    }
}

const vector<pair<double,double> >
IntonationHistogram::histogram(const vector<double> pitch)//, const vector<double> level)
{
    vector<pair<double,double> > out;
    const size_t numberOfBins   = m_freqs.size()-1;
    const double binWidth       = 1.0/m_nBPS;
    const double min            = m_freqs[0];
    const double hop            = binWidth/2.0;
    vector<double> count           = vector<double>(numberOfBins);
    
    // histo count
    for (size_t iPitch = 0; iPitch < pitch.size(); ++iPitch) {
        int bin = (int)((pitch[iPitch] - min) / binWidth);
        if (bin >= 0 && bin < numberOfBins) {
            count[bin] += 1; //level[iPitch];
        }
//        std::cerr << "pitch " << pitch[iPitch] << " " << binWidth << " " << bin << " "<< level[iPitch] << " "<< count[bin] << std::endl;
    }
    
    // histo area
    double histoArea = 0.0;
    for (size_t iBin = 0; iBin < numberOfBins; ++iBin) {
        histoArea += count[iBin]*binWidth;
    }
    
    for (size_t iBin = 0; iBin < numberOfBins; ++iBin) {
        pair<double,double> temp;
        temp.first  = m_freqs[iBin] + hop;
        temp.second = count[iBin]/histoArea;    // density
//        std::cerr << "count " << count[iBin] << " " << iBin << std::endl;

        out.push_back(temp);
    }
    
//    for (size_t iBin = 0; iBin < numberOfBins; ++iBin) {
//        std::cerr << "histo pair " << out[iBin].first << " "<< out[iBin].second << std::endl;
//    }
    
//    for (size_t iPitch = 0; iPitch < m_freqs.size(); ++iPitch) {
//        std::cerr << "m_freqs" << m_freqs[iPitch] << std::endl;
//    }
    
    return(out);
}

const std::vector<std::pair<double,double> >
IntonationHistogram::histogramFusion(vector<pair<double,double> > mpHisto, vector<pair<double,double> > npHisto, float fusionWeightStd)
{
    // create a weighting distribution
    boost::math::normal weightDistr = boost::math::normal(0, fusionWeightStd);
    
    // npHisto: note pitch track
    double tempMax = 0.0;
    double  tempMaxPitch = 0.0;
    for (size_t iBin = 0; iBin < npHisto.size(); ++iBin) {
        if (npHisto[iBin].second > tempMax) {
            tempMax = npHisto[iBin].second;
            tempMaxPitch = npHisto[iBin].first;
        }
    }
    
    // fusion
    vector<pair<double,double> > out;
    if (mpHisto.size() == npHisto.size()) {
        for (size_t iBin = 0; iBin < mpHisto.size(); ++iBin) {
            pair<double,double> temp;
            temp.first = mpHisto[iBin].first;
            double weightCoef = boost::math::pdf(weightDistr,std::fabs(temp.first - tempMaxPitch));
            temp.second = mpHisto[iBin].second * npHisto[iBin].second * weightCoef;
            out.push_back(temp);
        }
    }
    return(out);
}
