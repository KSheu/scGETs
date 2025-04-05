# scGETs Imputation
Single-cell Gene Expression Trajectory (scGETs) imputation package: a method for imputing stimulus-induced gene expression dynamics in single cells from time-series scRNASeq.

<img src="https://github.com/KSheu/scGETs/blob/main/GA_scGETs_protocol.png" width="350" height="350">

## Summary
Single-cell RNAseq (scRNAseq) measures cell-to-cell heterogeneous mRNA abundance but destroys the cell and precludes tracking of heterogeneous gene expression trajectories. Here we present an approach to impute single-cell gene expression trajectories (scGETs) from time-series scRNAseq measurements. We describe four main computational steps: dimensionality reduction, calculation of transition probability matrices, spline interpolation, and deconvolution to scGETs. Imputing scGETs can aid in studying heterogenous stimulus-responses over time, such as cancer cell responses to drugs, or immune cell responses to pathogens. 

For complete details on the use and execution of this protocol, please refer to Sheu et al, 2024. 


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
We can use the macrophage example data provided in the 'input' folder to run scREALTIME.
- Read in the Seurat object that contains annotated time-series scRNAseq data.
- Specify the timepoints to be used for the reconstruction, at least 4 timepoints. 
- Use the `getMetaData` function to obtain the metadata in the correct format. 
- Run the `scREALTIME` function and store the results. 

### Example Use
```
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC_scGETs_LPS.rds")
select_timepoints = c(0.0, 0.25, 1, 3, 8)
metadata = getMetaData(macro, stimulus = "LPS", timepoints= select_timepoints)
reconst = scREALTIME(input_obj = macro, metadata = metadata, timepoints = select_timepoints, stimulus = "LPS",
							num_archetypes = 20, num_sim_pts = 100, num_trajectories = 1000, 
							reduction = 'pca', consensus_measure = 'median', interpolant = 'spline', 
							data_norm = "ISnorm", prob_method = 'distance', distance_metric = 'euclidean' ,
							varFilter = T, exp_prob = 1) 
							
reconstructed_pc.traj = reconst$reconstructed_trajectories #Retrieve scGETs
```
The outputs of the imputation method are stored in a user-specified object, in this case `reconst`. scGETs are stored in the slot `reconst$reconstructed_trajectories`.

#### Plot the example output
```
reconstructed_pc.traj$stimulus = "LPS"
reconstructed_pc.traj$path_stim = paste0(reconstructed_pc.traj$path,"_", reconstructed_pc.traj$stimulus)
mat.numbers = reconstructed_pc.traj[,!grepl("time|path|stimulus|path_stim", colnames(reconstructed_pc.traj))]

# Rescale each gene column to range from 0-1
mat.numbers = apply(mat.numbers, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))) 

# Plot the result as lineplots
mat.numbers = cbind(mat.numbers, reconstructed_pc.traj[,grepl("time|path|stimulus|path_stim", colnames(reconstructed_pc.traj))])

gene = “Cxcl10”
ggplot(mat.numbers[grepl("",mat.numbers$stimulus),], aes(time,get(gene), group = path_stim)) +
  geom_vline(xintercept = c(0,0.25,1,3,8), linetype="dotted")+ 
  geom_line(aes(group = as.factor(path_stim)), alpha = 0.05)+
  theme_bw(base_size = 14)+theme(legend.position = "none")+ylab(gene)
```


### Input
Parameters specifying the input data are:
- `input_obj`: The Seurat Object containing the data.
- `metadata`: The table containing the metadata of the cells of interest.
- `timepoints`: A list of the timepoints in the dataset that are to be used for imputation
- `stimulus`: Name of the stimulus, in case the data object contains cells stimulated by multiple different stimuli. Impute scGETs for only one cell population stimulated with one stimulus at a time. 

### Parameters
Several parameters can be tuned by the user based on the characteristics of the dataset. 
- `num_archetypes` (default = 20): the number of cell subclusters at each timepoint. A smaller number will reduce the capture of outlier behavior in the resulting ensemble of trajectories but will speed up computational time. In general, a larger number of archetypes will better account for the single-cell distribution at each timepoint by increasing the search space of cell gene expression trajectories over time, at the expense of computational time and capturing unwanted technical noise. However, the number specified must be less than the smallest number of measured cells at any particular timepoint. 
- `num_sim_points` (default = 100): the number of interpolated data points across the timecourse. A larger number will provide a smoother trajectory. In specifying this number, consider the total length of the timecourse. The default value of 100 was used for an 8hr timecourse, equating a simulated point approximately every 5 minutes.
- `num_trajectories` (default = 1000): the number of total cells expected in the population. scREALTIME assumes that scRNAseq data is from time-series data over a short time course, with minimal cell death or division. Thus the total number of cells over time is assumed to be constant. Too small a number may not fully explore the space allowed by the transition probability matrix that links cell archetypes across timepoints. A recommendation is to specify a number close to the maximum number of cells measured over all the individual timepoints.

### Parameters (Advanced)
Several other hyperparameters may also be tuned for advanced customization of the protocol, listed below alongside their default values: 
- reduction = "pca": Dimensionality reduction to use ("pca","ica","nmf"). Must be an available reduction within the Seurat Object input. 
- consensus_measure = “median”: Method to define cell archetypes ("mean", "median").
- interpolant = “spline”: Method to interpolate timepoints ("spline", "linear").
- data_norm = "RNA": normalization from Seurat object to use (“RNA”, "SCT", "ISnorm").
- prob_method = "distance": Linkage probability based on a distance metric, or density of cells within archetypes (“distance”, “density”).
- distance_metric = "euclidean": Distance metric to use to identify cell archetype links over timepoints ("euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski").
- varFilter = T: Filters out zero variance genes prior to k-means clustering.
- exp_prob = 1: Raises the transition probability matrix to the power of exp_prob. Higher values make weak links weaker.


## Citation
Sheu, Katherine M., and Alexander Hoffmann. "Protocol: Imputation of stimulus-induced single-cell gene expression trajectories from time-series scRNAseq data." STAR Protocols (2025).

Sheu, Katherine M., Aditya Pimplaskar, and Alexander Hoffmann. “Single-Cell Stimulus-Response Gene Expression Trajectories Reveal the Stimulus Specificities of Dynamic Responses by Single Macrophages.” Molecular Cell 84, no. 21 (November 2024): 4095-4110.e6. https://doi.org/10.1016/j.molcel.2024.09.023.
