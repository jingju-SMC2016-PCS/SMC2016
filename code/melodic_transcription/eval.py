import os, sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'src'))

import subprocess
import evaluation
import vamp
import essentia.standard as es
import json, csv
import numpy as np

eval_noteTrans = evaluation.Evaluation()

# read ground truth note transcription
def readTranscription(monoNoteOut_filename):

        '''
        :param monoNoteOut_filename:
        :return: noteStartingTime, noteDurationTime, notePitch
        '''

        monoNoteOut         = np.loadtxt(monoNoteOut_filename,delimiter=',',usecols=[0,1,2])
        # print monoNoteOut
        noteStartingTime    = monoNoteOut[:,0]
        notePitch           = monoNoteOut[:,1]
        noteDurTime         = monoNoteOut[:,2]

        return [noteStartingTime, notePitch, noteDurTime]

# calculate note transcription with various parameters
def transcription_cal(audio_filename, monoNoteOut_filename):

    sr = 44100
    loader = es.MonoLoader(filename=audio_filename, downmix = 'mix', sampleRate = sr)
    audio = loader()

    # original pyin
    # data = vamp.collect(audio, sr, "pyin:pyin", output = 'notes')

    # modified pyin
    data = vamp.collect(audio, sr, "pyinbobigram:pyin-jingju", output = 'notes')

    # for note in data['list']:
    #     print note['duration'], note['timestamp'], note['values'][0]

    with open(monoNoteOut_filename, 'w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        for note in data['list']:
            csv_writer.writerow([note['timestamp'], note['values'][0], note['duration']])

# evaluate
def transcription_folder(folder):

    COnPOffList    = []
    COnPList       = []
    COnList        = []
    gt_lenList       = []
    gt_noteDurList   = []
    trans_lenList    = []

    filenames = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

    for filename in filenames:
        extension = os.path.splitext(filename)[1]
        if extension == '.wav':
            filenameRoot = os.path.splitext(filename)[0]
            filenameList = folder.split('/')
            audio_filename    = os.path.join(folder, filename)

            # estimate note transcription
            monoNoteOut_folder = './noteEstimatedOutput/'
            if not os.path.isdir(monoNoteOut_folder):
                os.makedirs(monoNoteOut_folder)

            monoNoteOut_filename = os.path.join(monoNoteOut_folder , filenameRoot+'.csv')

            print 'calculating audio file: ' + audio_filename

            transcription_cal(audio_filename, monoNoteOut_filename)

            print 'done! \n'

            # groundtruth
            gt_filename = os.path.join(folder, filenameRoot+'_melodicTrans.csv')
            gt = readTranscription(gt_filename)
            gt_noteDurList += gt[2].tolist()

            transcription = readTranscription(monoNoteOut_filename)

            # metrics
            COnPOff, COnP, COn = eval_noteTrans.coarseEval(gt, transcription)

            print 'number of COnPOff, COnP, COn'
            print COnPOff, COnP, COn
            print '\n'

            COnPOffList.append(COnPOff)
            COnPList.append(COnP)
            COnList.append(COn)

            gt_lenList.append(len(gt[0]))
            trans_lenList.append(len(transcription[0]))

            # print COnPOff, COnP, COn, len(gt[0]), len(transcription[0])

    return COnPOffList, COnPList, COnList, gt_lenList, gt_noteDurList, trans_lenList


folder_recordings                  = os.path.realpath('../../dataset/groundtruth/')

COnPOffList, COnPList, COnList, gt_lenList, gt_noteDurList, trans_lenList = transcription_folder(folder_recordings)
COnPOff, COnP, COn, gt_len, trans_len = sum(COnPOffList), sum(COnPList), sum(COnList), sum(gt_lenList), sum(trans_lenList)

# metrics
COnPOffF,COnPF,COnF, COnPOffP,COnPP,COnP, COnPOffR,COnPR,COnR = eval_noteTrans.metrics(COnPOff, COnP, COn, gt_len, trans_len)

print 'COnPOff, COnP, COn'
print 'F-measure'
print COnPOffF,COnPF,COnF 
print 'Precision'
print COnPOffP,COnPP,COnP 
print 'Recall'
print COnPOffR,COnPR,COnR
print 'groundtruth note number, mean duration, std duration'
print gt_len, np.mean(gt_noteDurList), np.std(gt_noteDurList)



