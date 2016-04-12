#Pitch contour segmentation
Segmentation of jingju pitch track into a chain of steady, transitory and vibrato segments.

##Usage
1. Make sure to clone the dataset folder which is required by this code.
2. Run optimization.py to obtain the complete parameter optimization results. There is a copy in "supplementary information" folder.
3. Run eval.py to evaluate the algorithm on test set.

##Dependencies
Essentia 2.0.1 http://essentia.upf.edu/  
Numpy 1.9.2  
Scipy 0.15.1  
Matplotlib 1.4.3  
Scikit-learn 0.16.1  
Statsmodels 0.6.1  
