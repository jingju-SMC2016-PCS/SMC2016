'''
interpolation of note transition distribution
'''

import numpy as np
import matplotlib.pyplot as plt

# distribution
intervals       = [-19, -15, -14, -12, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 17, 19]
intervalsFreq   = [4.80145964e-05,   1.92058386e-04,   1.44043789e-04,   2.40072982e-04,
   3.21697796e-03,   1.29639410e-03,   7.68233543e-03,   9.93902146e-03,
   3.36102175e-04,   5.13756182e-02,   1.76213569e-02,   1.08368944e-01,
   2.11456283e-01,   1.92538532e-02,   1.64257934e-01,   7.25020406e-03,
   1.86824795e-01,   1.18067893e-01,   1.99260575e-02,   4.20127719e-02,
   4.80145964e-05,   1.45004081e-02,   9.69894848e-03,   1.44043789e-03,
   3.55308014e-03,   9.60291929e-05,   8.16248139e-04,   9.60291929e-05,
   9.60291929e-05,   9.60291929e-05,   4.80145964e-05]

# interpolation
x               = np.linspace(-19,19,38*3+1)
y               = np.interp(x, intervals, intervalsFreq)

# normalisation constant
diffX           = x[1:] - x[:-1]
normConst       = 1/sum((y[:-1]+y[1:])*diffX/2.0)

yNorm           = y*normConst

# note transition distribution

print 'x', x
print 'y', y

# plt.stem(intervals,intervalsFreq)
plt.stem(x,y,'k')
plt.xlim(-5,5)

plt.xlabel("semitone distance from previous note")
plt.ylabel("probability")
plt.show()
