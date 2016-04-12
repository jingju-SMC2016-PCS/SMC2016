from math import *
import numpy as np
import matplotlib.pyplot as plt


class IntonationHistogram(object):

    def __init__(self):
        self.m_minFreq = 61.735
        self.m_nBPS = 5
        self.m_nPitch = 69 * self.m_nBPS  # 69 semi-tone, each semi-tone divided to 5, step is 20 cents
        self.m_freqs = np.zeros(self.m_nPitch, dtype=np.float64)
        for iPitch in range(self.m_nPitch):
            self.m_freqs[iPitch] = self.m_minFreq * pow(2, iPitch * 1.0 / (12 * self.m_nBPS))  # 0 to m_nPitch-1 positive pitch
            self.m_freqs[iPitch] = 12 * log(self.m_freqs[iPitch]/440.0)/log(2.0) + 69

    def histogram(self, pitch):
        '''
        Calcul the histogram of pitch array
        :param pitch:
        :return:
        '''

        hist, bin_edges = np.histogram(pitch, self.m_freqs, density=True)
        bin_edges = bin_edges[:-1] + 1/float(self.m_nBPS*2)

        # plt.figure()
        # plt.stem(bin_edges, hist)
        # plt.show()
        return hist, bin_edges
