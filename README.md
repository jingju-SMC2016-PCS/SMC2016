# SMC 2016 Pitch Contour Segmentation for Computer-aided Jingju Singing Training
This is the folder for containing the code, dataset and supplement information of SMC 2016 paper "Pitch Contour Segmentation for Computer-aided Jingju Singing Training"

## Code

The code used in this paper is in "code" folder, where you can find:

1. melodic_transcription: the modified pYIN algorithm code incorporated with the jingju bigram note transition probabilities + the binary for Mac OS X.
2. pitch_contour_segmentation: the python code for pitch contour segmentation without the "preliminary segmentation" step (which is already performed by using melodic_transcription).

## Dataset

The dataset used in this paper is in "dataset" folder. The a cappella singing audio recordings is not contained in this folder due to their large size, please contact the paper authors to request them (e-mail address will be given after the blind-reviewing process). In the folder you can find:

1. groundtruth: the ground truth annotation for 
 * melodic transcription (\*\_melodicTrans.csv, male\_12\_pos\_1 missing)
 * parameter optimization,
 * evaluating the StdCdLe thresholding (\*\_coarseSeg.csv) and the overall segmentation performance (\*\_refinedSeg.csv).


