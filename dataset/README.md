# SMC 2016 Pitch Contour Segmentation for Computer-aided Jingju Singing Training
This is the folder containing the dataset:

* file name: description (format)
* \*\_melodicTrans.csv: melodic transcription ground truth used for the evaluation (start_time pitch duration -).
* \*\_coarseSeg.csv: StdCdLe ground truth used for the parameter optimization and the evaluation (segmentation points).
* \*\_refinedSeg.csv: ground truth used for optimizing other parameters and the evaluation (start_time - duration).
* \*\_pitchtrack.csv: pitch track (contour) extracted by pYIN pitch-tracking algorithm (filename time pitch).
* \*\_monoNoteOut.csv: notes estimated by pYIN note-tracking algorithm (filename start_time duration pitch).
