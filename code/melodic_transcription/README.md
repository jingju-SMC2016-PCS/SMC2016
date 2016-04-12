#Melodic transcription for jingju singing
This folder contains the code for estimating the bigram note transition probabilities, and evaluating the transcription performance.

##Usage
* plotNoteTransitionProb.py: comparison plot of original pYIN note transition distribution and that estimated from jingju singing scores (dataset required).
* interpolateNoteTransitionProb.py: interpolation of bigram note transition probabilities for modified pYIN.
* eval.py: evaluation of the performance of melodic transcription performance (a cappella singing audio recordings required).

##Dependencies
* Essentia 2.0.1
* vamp python 1.1.0 https://pypi.python.org/pypi/vamp/1.0.0
* music21 http://web.mit.edu/music21/
* Numpy
* Scipy
* Matplotlib
