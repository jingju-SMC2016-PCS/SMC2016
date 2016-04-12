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
    datasetPath = os.path.realpath('../../dataset/')

    featureVecPath = './pYinOut/featureVec/'
    targetPath = './pYinOut/target/'

    # filenames for testing (evaluation)
    recordingNamesTest =    ['bcnRecording_001',
                            'bcnRecording_005',
                            'fem_01_neg_5',
                            'fem_01_pos_3',
                            'fem_10_pos_3',
                            'londonRecording_Dan-03',
                            'londonRecording_Laosheng-04',
                            'male_01_neg_3',
                            'male_01_pos_5',
                            'male_01_pos_6',
                            'male_13_pos_3']


    evaluation  = True                          # output evaluation results
    stdCdLe     = 0.22
    boundaryTh  = 0.02
    pitchTh     = 2.5                           #  5 semitones
    slopeTh     = 60.0                          #  contour combination slope difference threshold
    flatnoteTh  = 30.0                          #  threshold for deciding one note as flat pitch note


    
    ################################################## evaluation #########################################################
    # StdCdLe thresholding
    nc2 = nc.noteClass()
    nc2.noteSegmentationFeatureExtraction(datasetPath,
                                          featureVecPath,
                                          recordingNamesTest,
                                          segCoef=stdCdLe,
                                          predict=True)


    # pitch contour classification
    ttknn2 = ttknn.TrainTestKNN()
    ttknn2.predict(pitchContourClassificationModelName,
        featureVecPath,
        targetPath,
        recordingNamesTest)


    # refinement
    evalu1 = evalu.Evaluation()

    COnOffall, COnall, OBOnall, OBOffall,gt,st = 0,0,0,0,0,0        # evaluation metrics

    for rm in recordingNamesTest:

        #  filename declaration, 
        originalPitchtrackFilename              = os.path.join(datasetPath, rm+'_pitchtrack.csv')
        targetFilename                          = os.path.join(targetPath,  rm+'.json')
        refinedSegmentFeaturesFilename          = os.path.join(datasetPath, rm+'_refinedSegmentFeatures.json')
        representationFilename                  = os.path.join(datasetPath, rm+'_representation.json')
        figureFilename                          = os.path.join(datasetPath, rm+'_reprensentationContourFigure.png')
        regressionPitchtrackFilename            = os.path.join(datasetPath, rm+'_regression_pitchtrack.csv')
        refinedSegmentationGroundtruthFilename  = os.path.join(datasetPath, rm+'_refinedSeg.csv')

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

        COnOffall   += rsm1.evaluationMetrics[0]
        COnall      += rsm1.evaluationMetrics[1]
        OBOnall     += rsm1.evaluationMetrics[2]
        OBOffall    += rsm1.evaluationMetrics[3]
        gt          += rsm1.evaluationMetrics[4]
        st          += rsm1.evaluationMetrics[5]

    COnOffF,COnF,COnOffP,COnP,COnOffR,COnR,OBOnRateGT,OBOffRateGT = evalu1.metrics(COnOffall,COnall,OBOnall,OBOffall,gt,st)

    print 'Parameters: Th_time, Th_pitch, Th_slope, Th_flat pitch'
    print boundaryTh,pitchTh,slopeTh,flatnoteTh
    print 'Metrics of COnOff: F-measure, Precision, Recall'
    print COnOffF,COnOffP,COnOffR
    print 'Metrics of COn: F-measure, Precision, Recall'
    print COnF,COnP,COnR
    

    print("--- %s seconds ---" % (time.time() - overall_start_time))
