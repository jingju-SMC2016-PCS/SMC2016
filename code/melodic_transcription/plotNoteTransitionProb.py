'''
calculate note distribution from score
'''

import os, sys
from os import listdir
from os.path import isfile, join, splitext
import matplotlib.pyplot as plt
import numpy as np
import collections
from scipy.stats import norm

sys.path.append(join(os.path.dirname(os.path.realpath(__file__)),'src'))

import musicXMLparser

##########################################
#  jingju note transition probabilities  #
##########################################

scoreFolder = os.path.realpath('../../dataset/scores/')
onlyfiles = [f for f in listdir(scoreFolder) if isfile(join(scoreFolder, f))]

scoreTotalNum       = 0
labelFontsize 		= 18

intervals_total = []
for filename in onlyfiles:
    filenameRoot, fileExtension = splitext(filename)
    if fileExtension == '.xml':
        # get intervals of the singing track
        intervals = musicXMLparser.getIntervals(join(scoreFolder, filename))
        intervals_total += intervals
        scoreTotalNum += 1

# find unique elements
intervalsSet    = set(intervals_total)
intervalsUnique = list(intervalsSet)

intervalsFreq = {}
# construct the dictionary
for iu in intervalsUnique:
    intervalsFreq[iu] = 0

# fill the dictionary
for iu in intervalsUnique:
    for ii in intervals_total:
        if iu == ii:
            intervalsFreq[iu] += 1

# sort the dictionary by keys
oIntervalsFreq = collections.OrderedDict(sorted(intervalsFreq.items()))

# keys and values
interFreqKeys   = oIntervalsFreq.keys()
interFreqValues = oIntervalsFreq.values()

# fit intervals to norm distribution
interMean, interStd = norm.fit(intervals_total)
interX              = range(min(interFreqKeys),max(interFreqKeys)+1)
interPdf            = norm.pdf(interX,interMean,interStd)
interFreqValuesNorm = np.array(interFreqValues)/float(sum(interFreqValues))

# print oIntervalsFreq
# print 'interFreqKeysNorm', interFreqKeys                # intervals
# print 'interFreqValuesNorm', interFreqValuesNorm        # intervals frequency
# print interPdf
# print 'score total num: ', scoreTotalNum
# print 'intervals total num: ', len(intervals_total)


########################################
#  pYIN note transition probabilities  #
########################################

sigma2Note = 0.7  # note transition probability standard deviation used in pYin

noteDistanceDistr = norm(loc=0, scale=sigma2Note)

midIndex = int(21/2) # x middle index

x = np.linspace(-10,10,21)

y = norm.pdf(x)

y[midIndex - 1] = 0
y[midIndex + 1] = 0

##########
#  plot  #
##########

# subplot pyin note transition probabilities
plt.subplot(2, 1, 1)
plt.stem(x,y,'k')
plt.axis([-10,10, 0, max(y)*1.2])
plt.ylabel("probability", fontsize = labelFontsize)

# subplot jingju note transition probabilities
plt.subplot(2, 1, 2)
plt.stem(interFreqKeys,interFreqValuesNorm,'k')
plt.xlim(-10,10)
plt.ylabel("probability", fontsize = labelFontsize)
plt.xlabel("semitone distance from previous note", fontsize = labelFontsize)
plt.show()



