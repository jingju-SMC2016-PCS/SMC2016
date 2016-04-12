'''
This is the code for evaluation.

the coarse evaluation: groundtruth is corrected on the pYIN algorithm by adding stable part segmentation.
the fine evaluation: groundtruth is the vibrato, flat, ascending and descending segmentation which is not based on pYIN segmentation.

'''

import numpy as np


class Evaluation(object):

    def precision(self,c,groundtruth):
        return c/float(groundtruth)

    def recall(self,c,seg):
        return c/float(seg)

    def Fmeasure(self,p,r):
        return 2*(p*r)/(p+r)

    def errorRateGT(self,e,groundtruth):
        return e/float(groundtruth)

    def errorRateTR(self,e,seg):
        return e/float(seg)

    def errorRatio(self,rateGT,rateTR):
        if not rateGT:
            return None
        else:
            return rateTR/rateGT

    def offsetTh(self,groundtruth):

        thOffset = []

        for gt in groundtruth[2]:
            twentyP = gt*0.2
            thOffset.append(max(twentyP,0.05))

        return thOffset

    def coarseEval(self, groundtruth, seg):
        '''
        :param groundtruth: a list [startTime, pitch, endTime,...]
        :param seg: segementation, a list [startTime, pitch, endTime,...]
        :return: COnOff, COn, OBOn, OBOff
        '''

        COnPOff = 0             # correct onset offset pitch
        COnP = 0                # correct onset pitch
        COn = 0                 # correct onset

        thOnset = 0.05          # 50 ms
        thOffset = self.offsetTh(groundtruth)

        startTime_s = np.array(seg)[0]
        endTime_s   = np.array(seg)[0] + np.array(seg)[2]
        pitch_s     = 12 * np.log(np.array(seg)[1]/440)/np.log(2.) + 69

        startTime_gt = np.array(groundtruth)[0]
        endTime_gt   = np.array(groundtruth)[0] + np.array(groundtruth)[2]
        pitch_gt     = 12 * np.log(np.array(groundtruth)[1]/440)/np.log(2.) + 69

        # print 'startTime_s',    startTime_s
        # print 'endTime_s',      endTime_s
        # print 'pitch_s',        pitch_s
        #
        # print 'startTime_gt',   startTime_gt
        # print 'endTime_gt',     endTime_gt
        # print 'pitch_gt',       pitch_gt

        for si in range(len(seg[0])):
            for gti in range(len(groundtruth[0])):
                # COn
                if abs(startTime_gt[gti] - startTime_s[si]) < thOnset:
                    COn += 1
                    # COnOff
                    if abs(pitch_gt[gti] - pitch_s[si]) < 0.5:
                        COnP += 1
                        if abs(endTime_gt[gti] - endTime_s[si]) < thOffset[gti]:
                            COnPOff += 1
                    break       # break the first loop if s have been already mapped to


        return COnPOff, COnP, COn
        #print groundtruth, seg


    def metrics(self,COnPOff, COnP, COn, groundtruthLen,segLen):

        COnPOffP = self.precision(COnPOff,groundtruthLen)
        COnPOffR = self.recall(COnPOff,segLen)
        COnPOffF = self.Fmeasure(COnPOffP,COnPOffR)

        COnPP = self.precision(COnP,groundtruthLen)
        COnPR = self.recall(COnP,segLen)
        COnPF = self.Fmeasure(COnPP,COnPR)

        COnP_ = self.precision(COn,groundtruthLen)
        COnR_ = self.recall(COn,segLen)
        COnF_ = self.Fmeasure(COnP_,COnR_)

        return COnPOffF,COnPF,COnF_, COnPOffP,COnPP,COnP_, COnPOffR,COnPR,COnR_




