# -*- coding: utf-8 -*-

import sys, os
import matplotlib.pyplot as plt
import time
import subprocess
import json
import shutil
import random

# add src path
dir = os.path.dirname(os.path.realpath(__file__))
srcpath = dir+'/src'
sys.path.append(srcpath)

import pitchtrackSegByNotes
import essentia.standard as ess
import numpy as np
import noteClass as nc
import featureVecTarget as fvt
import trainTestKNN as ttknn
import refinedSegmentsManipulation as rsm
import evaluation as evalu
from vibrato import vibrato

if __name__ == "__main__":

    overall_start_time = time.time()  # starting time

    # classification model path
    pitchContourClassificationModelName = './pYinOut/train/model/pitchContourClassificationModel.pkl'

    # ground truth path
    datasetPath         = os.path.realpath('../../dataset/groundtruth/')

    featureVecPath      = './pYinOut/featureVec/'
    targetPath          = './pYinOut/target/'

    evaluation          = True

    # filenames for parameter optimization (training dataset)
    recordingNamesTraining =     ['bcnRecording_003',
                                'bcnRecording_004',
                                'bcnRecording_007',
                                'bcnRecording_008',
                                'fem_01_neg_1',
                                'fem_01_neg_3',
                                'fem_01_pos_1',
                                'fem_01_pos_5',
                                'fem_01_pos_7',
                                'fem_07_pos_1',
                                'fem_10_pos_1',
                                'fem_11_pos_1',
                                'londonRecording_Dan-01',
                                'londonRecording_Dan-02',
                                'londonRecording_Dan-04',
                                'londonRecording_Laosheng-01',
                                'londonRecording_Laosheng-02',
                                'londonRecording_Laosheng-03',
                                'male_01_neg_1',
                                'male_01_neg_2',
                                'male_01_neg_4',
                                'male_01_neg_5',
                                'male_01_pos_1',
                                'male_01_pos_2',
                                'male_01_pos_3',
                                'male_01_pos_4',
                                'male_02_neg_1',
                                'male_12_neg_1',
                                'male_12_pos_1',
                                'male_13_pos_1'
                                ]

    
    ################################################ parameter optimization ###################################################
    # uncomment it only when needs evaluation

    # StdCdLe grid search
    nc2 = nc.noteClass()
    with open('./optimization_results/stdCdLeGridSearch.txt', "w") as outfile:
        for sc in np.linspace(0.01,1.0,100):
            COnOffF,COnF, COnOffP,COnP,COnOffR, COnR,OBOnRateGT,OBOffRateGT = \
                nc2.noteSegmentationFeatureExtraction(datasetPath,
                                                    featureVecPath,
                                                    recordingNamesTraining,
                                                    segCoef=sc,
                                                    predict=True,
                                                    evaluation=True)
            outfile.write(str(sc)+'\t'+str([COnOffF,COnF, COnOffP,COnP,COnOffR, COnR,OBOnRateGT,OBOffRateGT])+'\n')

    # pitch contour classification
    ttknn2 = ttknn.TrainTestKNN()
    ttknn2.predict(pitchContourClassificationModelName,
        featureVecPath,
        targetPath,
        recordingNamesTest)
 
    # other parameter grid search

    evalu1 = evalu.Evaluation()
    with open('./optimization_results/otherParametersGridSearch.txt', "w") as outfile:
        #  grid search
        for boundaryTh in [0.01,0.02,0.05,0.1,0.2,0.5]:
            for pitchTh in [0.1,0.2,0.5,1,2,5,10]:
                pitchTh = pitchTh/2.0   # semitone
                for slopeTh in range(10,100,10):
                    for flatnoteTh in range(10,100,10):

                        start_time = time.time()  # starting time

                        if evaluation:
                            COnOffall, COnall, OBOnall, OBOffall,gt,st = 0,0,0,0,0,0        # evaluation metrics

                        for rm in recordingNamesTraining:

                            print 'processing file: ', rm, ' boundaryTh, ', boundaryTh, ' pitchTh,', pitchTh, ' slopeTh, ', slopeTh, ' flatenoteTh, ', flatnoteTh
                            
                            #  filename declaration
                            originalPitchtrackFilename              = os.path.join(datasetPath, rm+'_pitchtrack.csv')
                            targetFilename                          = os.path.join(targetPath,  rm+'.json')
                            refinedSegmentFeaturesFilename          = os.path.join(datasetPath, rm+'_refinedSegmentFeatures.json')
                            representationFilename                  = os.path.join(datasetPath, rm+'_representation.json')
                            figureFilename                          = os.path.join(datasetPath, rm+'_reprensentationContourFigure.png')
                            regressionPitchtrackFilename            = os.path.join(datasetPath, rm+'_regression_pitchtrack.csv')

                            if evaluation:
                                refinedSegmentationGroundtruthFilename  = os.path.join(datasetPath, rm+'_refinedSeg.csv')
                            else:
                                refinedSegmentationGroundtruthFilename = None

                            rsm1 = rsm.RefinedSegmentsManipulation()
                            rsm1.process(refinedSegmentFeaturesName     = refinedSegmentFeaturesFilename,
                                         targetFilename                 = targetFilename,
                                         representationFilename         = representationFilename,
                                         figureFilename                 = figureFilename,
                                         regressionPitchtrackFilename   = regressionPitchtrackFilename,
                                         originalPitchtrackFilename     = originalPitchtrackFilename,
                                         refinedSegGroundtruthFilename  = refinedSegmentationGroundtruthFilename,
                                         evaluation                     = evaluation,
                                         slopeTh                        = slopeTh,
                                         flatnoteTh                     = flatnoteTh,
                                         boundaryTh                     = boundaryTh,
                                         pitchTh                        = pitchTh)

                            #  evaluation metrics collection
                            if evaluation:                                          
                                COnOffall   += rsm1.evaluationMetrics[0]
                                COnall      += rsm1.evaluationMetrics[1]
                                OBOnall     += rsm1.evaluationMetrics[2]
                                OBOffall    += rsm1.evaluationMetrics[3]
                                gt          += rsm1.evaluationMetrics[4]
                                st          += rsm1.evaluationMetrics[5]

                        if evaluation:
                            COnOffF,COnF,COnOffP,COnP,COnOffR,COnR,OBOnRateGT,OBOffRateGT = evalu1.metrics(COnOffall,COnall,OBOnall,OBOffall,gt,st)
                            outfile.write(str(boundaryTh)+'\t'+str(pitchTh)+'\t'+
                                          str(slopeTh)+'\t'+str(flatnoteTh)+'\t'+
                                          str(COnOffF)+'\t'+str(COnF)+'\t'+
                                          str(COnOffP)+'\t'+str(COnP)+'\t'+
                                          str(COnOffR)+'\t'+str(COnR)+'\t'+
                                          str(OBOnRateGT)+'\t'+str(OBOffRateGT)+'\n')
                            
                            print 'Parameters: Th_time, Th_pitch, Th_slope, Th_flat pitch'
                            print boundaryTh,pitchTh,slopeTh,flatnoteTh
                            print 'Metrics of COnOff: F-measure, Precision, Recall'
                            print COnOffF,COnOffP,COnOffR
                            print 'Metrics of COn: F-measure, Precision, Recall'
                            print COnF,COnP,COnR

                        print("--- %s seconds ---" % (time.time() - start_time))
    
    

    print("--- %s seconds ---" % (time.time() - overall_start_time))
