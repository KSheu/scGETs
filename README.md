# scGETs Imputation
Single-cell Gene Expression Trajectory (scGETs) imputation package: a method for imputing stimulus-induced gene expression dynamics in single cells from time-series scRNASeq.

<img src="https://github.com/KSheu/scGETs/blob/main/GA_scGETs_protocol.png" width="350" height="350">


## Dependencies
- Seurat
- Nbclust
- matrixStats
- gridExtra

## Install
Tested compatilbility with R version >=4.2.1.\
Install with: 
```
library(devtools)
install_github("ksheu/scGETs")
library(scGETs)
```

## Quick Start
We can use the macrophage example data provided in the 'output' folder to run scREALTIME.
- Read in the Seurat object that contains annotated time-series scRNAseq data.
- Specify the timepoints to be used for the reconstruction, at least 4 timepoints. 
- Use the getMetaData function to obtain the metadata in the correct format. 
- Run the scREALTIME function and store the results. 

### Example Use
```
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
select_timepoints = c(0.0, 0.25, 1, 3, 8)
metadata = getMetaData(macro, stimulus = "LPS", timepoints= select_timepoints)
reconst = scREALTIME(input_obj = macro, metadata = metadata, timepoints = select_timepoints, stimulus = "LPS",
							num_archetypes = 20, num_sim_pts = 100, num_trajectories = 1000, 
							reduction = 'pca', consensus_measure = 'median', interpolant = 'spline', 
							data = "ISnorm", prob_method = 'distance', distance_metric = 'euclidean' ,
							varFilter = T, exp_prob = 1) 
							

```
The outputs of the imputation method are stored in a user-specified object, in this case “reconst”. scGETs are stored in the slot `reconst$reconstructed_trajectories`.

### Input
Parameters specifying the input data are:
- input_obj: The Seurat Object containing the data.
- metadata: The table containing the metadata of the cells of interest.
- timepoints: A list of the timepoints in the dataset that are to be used for imputation
- stimulus: Name of the stimulus, in case the data object contains cells stimulated by multiple different stimuli. Impute scGETs for only one cell population stimulated with one stimulus at a time. 


### Parameters
Several parameters can be tuned by the user based on the characteristics of the dataset. 
- num_archetypes (default = 20): The number of cell subclusters at each timepoint. A larger number will better account for the single-cell distribution at each timepoint. A smaller number will reduce the capture of outlier behavior in the resulting ensemble of trajectories but will speed up computational time. Number specified must be less than the smallest number of measured cells at any particular timepoint.
- num_sim_points (default = 100): The number of interpolated data points across the timecourse. A larger number will provide a smoother trajectory.
- num_trajectories (default = 1000): The number of total cells expected in the population. scREALTIME assumes that scRNAseq data is from time-series data over a short time course, with minimal cell death or division. Thus the total number of cells over time is assumed to be constant. Too small a number may not fully explore the space allowed by the transition probability matrix that links cell archetypes across timepoints.

### Parameters (Advanced)
Several other hyperparameters may also be tuned for advanced customization of the protocol, listed below alongside their default values: 
- reduction = "pca": Dimensionality reduction to use ("pca","ica","nmf").
- consensus_measure = “median”: Method to define cell archetypes ("mean", "median").
- interpolant = “spline”: Method to interpolate timepoints ("spline", "linear").
- data = "RNA": normalization from Seurat object to use (“RNA”, "SCT", "ISnorm").
- prob_method = "distance": Linkage probability based on a distance metric, or density of cells within archetypes (“distance”, “density”).
- distance_metric = "euclidean": Distance metric to use to identify cell archetype links over timepoints. Must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
- varFilter = T: Filters out zero variance genes prior to k-means clustering.
- exp_prob = 1: Raises the transition probability matrix to the power of exp_prob. Higher values make weak links weaker.

