# Projection-Clustering
This repository contains supporting codes for manuscript *Bayesian clustering using random effects models and predictive projections*.
  
## Examples
1. Synthetic dataset 
    1. Data observations![Example 1 data](eg1_synthetic/plotSeriesEg1.png)
    2. Clustering probabilities ![Example 1 clustering](eg1_synthetic/plotClusterProbEg1.png)
2. [Crop image](https://www.cs.ucr.edu/%7Eeamonn/time_series_data_2018/)
4. [DNA synchrony of yeast cells](http://genome-www.stanford.edu/cellcycle/)
5. [EEG signals during sleep](https://physionet.org/content/capslpdb/1.0.0/)
6. [Activity recognition from accelerometer data](https://archive.ics.uci.edu/ml/datasets/Activity+Recognition+from+Single+Chest-Mounted+Accelerometer)

## Note on benchmark methods
All benchmarks make use of existing R packages, except for the variational computation method adapted from paper *Variational Approximation for Mixtures of Linear Mixed Models* by Tan & Nott,  of which codes are available on [supplemental materials](https://doi.org/10.1080/10618600.2012.761138).

## Note on application
To recreate analysis described in manuscript, enter example-specific subfolder and run code.
