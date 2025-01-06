# scGETs Imputation
Single-cell Gene Expression Trajectory (scGETs) imputation package: a method for imputing stimulus-induced gene expression dynamics in single cells from time-series scRNASeq.

<img src="https://github.com/KSheu/scGETs/blob/main/GA_scGETs_protocol.png" width="100" height="100">


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
- Run the scGETs function and store the results. 

### Example Use
```
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
select_timepoints = c(0.0, 0.25, 1, 3, 8)
metadata = getMetaData(macro, stimulus = "LPS", timepoints= select_timepoints)
reconst = scGETs(macro = macro, metadata = metadata, timepoints = select_timepoints,
							num_archetypes = 20, data = "ISNorm",
							num_trajectories = 1000, num_sim_pts = 100,
							reduction = 'pca', stimulus = "LPS", consensus_measure = 'median',
							interpolant = 'spline', prob_method = 'distance', distance_metric = 'euclidean' ,
							varFilter = T, exp_prob = 1) 
							

```
scREALTIME trajectories are stored in the slot `reconst$reconstructed_trajectories`.

### Parameters
Several parameters can be tuned by the user based on the characteristics of the dataset. 
- num_archetypes: controls the number of cell subclusters at each timepoints. A larger number will better account for the single-cell distribution at each timepoint. A smaller number will minimize the capture of outlier behavior in the resulting ensemble of trajectories. 
- num_sim_points: number of interpolated data points across the timecourse. A larger number will provide a smoother trajectory. 
- num_trajectories: number of total cells expected in the population. scREALTIME assumes that scRNAseq data is from a short timecourse, with minimal cell death or division. Too small a number may not fully explore the space allowed by the transition probability matrix that links cell archetypes across timepoints.
