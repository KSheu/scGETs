# single-cell response trajectory imputation algorithm
# auxilliary functions for plotting and data extraction
# contact: katherinesheu[at]ucla[dot]edu


#' scREALTIME - single-cell Reconstruction of Expression using Archetype Linkage of TIme-series MEasurements
#' Input a Seurat object of single cell data and returns list object containing single cell trajectories
#' @param macro Seurat object of single cell data with multiple measured timepoints
#' @param metadata Metadata frame of the measured single cells, obtain with getMetaData
#' @param timepoints Timepoints to use from the measured data
#' @param stimulus Specify one stimulus to impute if desired
#' @param num_archetypes = 20 Number of cell archetypes expected per timepoint (eg. 10, 20)
#' @param num_sim_pts = 100 Number of simulated timepoints to output
#' @param num_trajectories = 1000 Number of simulated trajectories expected
#' @param reduction = "pca" Dimensionality reduction to use, must be available in Seurat Object ("pca","ica","nmf")
#' @param consensus_measure Method to define cell archetypes ("mean", "median")
#' @param interpolant Method to interpolate timepoints ("spline","linear")
#' @param data = "RNA" normalization from Seurat object to use ("RNA", "SCT", "ISnorm")
#' @param prob_method = "distance" Linkage probability based on a 'distance' metric, or 'density' of cells in archetypes("distance", "density")
#' @param distance_metric = "euclidean" Distance metric to use to identify cell archetype links over timepoints. Must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param varFilter = T Filters out zero variance genes prior to kmeans clustering.
#' @param exp_prob = 1 Raises the transition probability matrix to the power of exp_prob. Higher values make weak links weaker.
#' @return trajectoryObject

# scREALTIME
# Input a Seurat object of single cell data and returns list object containing single cell trajectories
scREALTIME = function(input_obj, metadata, num_archetypes=20, timepoints, num_trajectories = 1000, num_sim_pts = 100,
                      reduction = "pca", stimulus, consensus_measure = "median", interpolant="spline", data = "RNA", prob_method = 'distance', distance_metric = 'euclidean', varFilter = T, exp_prob = 1){

  # Retrieve desired data subset for each timepoint
  cells_by_timept <- list()
  for(i in timepoints){
    index = paste("time_",i, "hr", sep = "")
    cells_by_timept[[index]] <- rownames(metadata)[metadata$timept_num == i]
  }

  # Retrieve RNA data based on specified normalization
  if(data == 'RNA'){
    RNA <- as.data.frame(input_obj@assays$RNA@data)
  }
  else if(data == 'ISnorm'){
    RNA <- as.data.frame(input_obj@assays$ISnorm@data)
  }
  RNA <- t(RNA)

  # Filter RNA data to only include cells in the metadata
  RNA <- RNA[rownames(RNA) %in% rownames(metadata),]

  # Perform k-means clustering for each timepoint
  clusterings <- list()
  zero_var_inds <- list()
  for(i in timepoints){
    index = paste("time_", i, "hr", sep = "")
    data <- RNA[rownames(RNA) %in% cells_by_timept[[index]],]
    zero_var_inds[[index]] <- resample::colVars(data) == 0

    if(varFilter){
      clusterings[[index]] <- kmeans(data[,!zero_var_inds[[index]]], centers = num_archetypes, iter.max = 50)
    }else{
      clusterings[[index]] <- kmeans(data, centers = num_archetypes, iter.max = 50)
    }

    # Track cluster counts for each timepoint
    if(i == timepoints[1]){
      cluster_counts <- as.data.frame(table(clusterings[[index]]$cluster))
    }else{
      cluster_counts <- cbind(cluster_counts, table(clusterings[[index]]$cluster))
    }
  }

  # Clean up cluster counts data frame
  col_omit <- c(1:length(timepoints))*2 - 1
  cluster_counts <- cluster_counts[, -col_omit]
  rownames(cluster_counts) <- c(1:num_archetypes)
  rownames(cluster_counts) <- paste("bin", rownames(cluster_counts), "")
  colnames(cluster_counts) <- timepoints
  colnames(cluster_counts) <- paste("time_", colnames(cluster_counts), sep="")

  # Calculate and store cluster density from cluster_counts table
  cluster_densities <- as.data.frame(prop.table(as.matrix(cluster_counts), 2))
  cluster_densities <- cluster_densities^(exp_prob)

  # Retrieve PCA (or other dimensionality reduction) results from Seurat object
  if(toupper(reduction) == 'PCA'){
    pcscores = input_obj[['pca']]@cell.embeddings
  }else if(toupper(reduction) == 'NMF'){
    pcscores = input_obj[['NMF']]@cell.embeddings
  }else if(toupper(reduction) == 'ICA'){
    pcscores = input_obj[['ica']]@cell.embeddings
  }


  # Identify cell archetypes in PC space
  cell_cluster_df <- matrix(nrow = 0, ncol = 2)
  for(i in timepoints){
    index = paste("time_", i, "hr", sep = "")
    cell_cluster_df <- c(cell_cluster_df, clusterings[[index]]$cluster)
  }
  cell_cluster_df <- as.data.frame(cell_cluster_df)

  pcscores <- as.data.frame(pcscores)
  pcscores_stim <- pcscores[rownames(pcscores) %in% rownames(metadata),]
  pcscores_stim$timept <- metadata$timept_num[match(rownames(pcscores_stim), rownames(metadata))]
  pcscores_stim$bin <- cell_cluster_df$cell_cluster_df[match(rownames(pcscores_stim), rownames(cell_cluster_df))]
  pcscores_stim$timebin_tag <- paste(pcscores_stim$timept, pcscores_stim$bin, sep = "_")

  # Calculate mean/median PC scores for each timepoint-bin combination
  if(consensus_measure == 'mean'){
    aggreg_pcscores <- aggregate(pcscores_stim[, 1:(ncol(pcscores_stim)-3)], list(pcscores_stim$timebin_tag), mean)
  }else if(consensus_measure == 'median'){
    aggreg_pcscores <- aggregate(pcscores_stim[, 1:(ncol(pcscores_stim)-3)], list(pcscores_stim$timebin_tag), median)
  }
  aggreg_pcscores$timept <- sapply(strsplit(aggreg_pcscores$Group.1, split = "_"), `[`, 1)
  aggreg_pcscores$bin <- sapply(strsplit(aggreg_pcscores$Group.1, split = "_"), `[`, 2)
  col_orders = c(1,(ncol(aggreg_pcscores)), (ncol(aggreg_pcscores)-1), 2:(ncol(aggreg_pcscores)-2))
  aggreg_pcscores <- aggreg_pcscores[order(aggreg_pcscores$timept),col_orders]
  aggreg_pcscores <- aggreg_pcscores[order(as.numeric(aggreg_pcscores$bin)),]
  rownames(aggreg_pcscores) = aggreg_pcscores$Group.1

  # Generate transition probability matrices based on distance between clusters
  distances = as.matrix(dist(aggreg_pcscores[,4:ncol(aggreg_pcscores)],method = distance_metric))
  rownames(distances) = aggreg_pcscores$Group.1
  colnames(distances) = aggreg_pcscores$Group.1

  # Generate random walks based on distance between clusters
  walks <- list()
  walk_probs <- matrix(nrow = num_trajectories, ncol = 2)
  walk_probs[,1] = c(1:num_trajectories)

  if(prob_method == 'distance'){
    for(i in 1:num_trajectories){
      path <- c()
      for(j in 1:length(timepoints)){
        if(j==1){
          next_value = sample(c(1:num_archetypes), size = 1, prob = cluster_densities[,j])
        }
        else{
          prev = path[length(path)]
          row = (num_archetypes)*(j-2) + prev
          cols = (num_archetypes*(j-1) + 1):(num_archetypes*j)
          probs = distances[row, cols]
          probs = (1/probs)^exp_prob
          next_value = sample(c(1:num_archetypes), size = 1, prob = probs)
        }
        path <- c(path, next_value)
      }
      walks[[i]] <- path
    }
  }else if(prob_method == 'density'){
    for(i in 1:num_trajectories){
      cur_prob = 1
      path <- c()
      for(j in 1:length(timepoints)){
        next_value = sample(c(1:num_archetypes), size = 1, prob = cluster_densities[,j])
        cur_prob = cur_prob * cluster_densities[next_value,j]
        path <- c(path, next_value)
      }
      walks[[i]] <- path
      walk_probs[i,2] <- cur_prob
    }
  }else{
    print('Invalid method for random walk probabilities')
  }

  # Check number of unique walks
  num_unique_traj = length(unique(sapply( walks, paste0, collapse="")))
  print(paste(length(unique(sapply( walks, paste0, collapse=""))), "unique walks/trajectories"))

  # Make table of unique walks
  walk_frequencies <- as.data.frame(table(sapply( walks, paste0, collapse="")))
  colnames(walk_frequencies)[1] = 'Walk'

  # Spline fitting across the random walks
  # Generate simulated time points and manually add the measured timepoints to num_sim_pts
  spline_pts <- list()
  for(pc in 1:(ncol(aggreg_pcscores)-3)){
    id = paste("pc", pc, sep="_")
    spline_pts[[id]] <- matrix(nrow = num_trajectories, ncol = length(timepoints))
    for(tr in 1:num_trajectories){
      path = walks[[tr]]
      for(j in 1:length(timepoints)){
        time = timepoints[j]
        bin = path[j]
        spline_pts[[id]][tr, j] <- aggreg_pcscores[round(as.numeric(aggreg_pcscores$timept),3) == round(time,3) & aggreg_pcscores$bin == bin, pc+3]
      }
    }
  }

  # Handle duplicated rows in spline points
  dup_rows = duplicated(spline_pts[[1]])
  walk_probs = walk_probs[!dup_rows, ]

  for(i in 1:length(spline_pts)){
    spline_pts[[i]] <- spline_pts[[i]][!dup_rows,]
  }

  # Generate simulated time points
  sim_times = seq(min(timepoints),max(timepoints),length.out = num_sim_pts)
  for(i in timepoints){
    if(!(i %in% sim_times)){
      sim_times <- c(sim_times, i)
      num_sim_pts = num_sim_pts + 1
    }
  }
  sim_times <- sort(sim_times)

  simulated = matrix(ncol = (ncol(aggreg_pcscores)-3), nrow = num_sim_pts*dim(walk_frequencies)[1])
  simulation_mat_time = rep(sim_times, dim(walk_frequencies)[1])

  # Perform spline interpolation using an interpolation method
  for(i in 1:(ncol(aggreg_pcscores)-3)){
    preds = c()
    for(walk in 1:nrow(spline_pts[[i]])){
      if(interpolant == 'spline'){
        spline <- smooth.spline(timepoints, spline_pts[[i]][walk,], cv = T, all.knots = T)
        preds <- c(preds, predict(spline, sim_times)$y)
      } else if(interpolant == 'linear'){
        fun = approxfun(timepoints, spline_pts[[i]][walk,], rule = 2)
        preds = c(preds, fun(sim_times))
      } else if(interpolant == 'loess'){
        loe = loess(spline_pts[[i]][walk,] ~ timepoints, span = 0.50, se = F)
        preds = c(preds, predict(loe, sim_times))
      }else {
        print('Unspecified interpolant')
      }
    }
    simulated[,i] = preds
  }

  # Cast back to gene expression space
  loadings <- input_obj[[reduction]]@feature.loadings

  if(toupper(reduction) == 'NMF'){
    loadings = input_obj[['NMF']]@feature.loadings
  }
  if(toupper(reduction) == 'ICA'){
    loadings = input_obj[['ica']]@feature.loadings
  }
  reconstructed_pc <- t(loadings %*% t(simulated))
  reconstructed_pc <- as.data.frame(reconstructed_pc)
  reconstructed_pc$time <- simulation_mat_time
  reconstructed_pc$path <- rep(1:dim(walk_frequencies)[1], each = num_sim_pts)

  # Return result
  toRet = list()
  toRet[['reconstructed_trajectories']] = reconstructed_pc
  toRet[['cluster_densities']] = cluster_densities
  toRet[['metadata']] = metadata
  toRet[['number_unique_trajectories']] = num_unique_traj
  toRet[['probability_of_walks']] = walk_probs

  call = list()
  call[['Seurat_obj']] = input_obj
  call[['metadata_name']] = metadata
  call[['num_archetypes']] = num_archetypes
  call[['timepoints']] = timepoints
  call[['num_trajectories']] = num_trajectories
  call[['num_sim_pts']] = num_sim_pts
  call[['reduction']] = reduction
  call[['stimulus']] = stimulus
  call[['consensus_measure']] = consensus_measure
  toRet[['call']] = call

  return(toRet)
}

#' Input a Seurat object of single cell data and returns list object containing the metadata
#' @param macro Seurat object of single cell data with multiple measured timepoints
#' @param stimulus Name of the stimulus wanted, assumes metadata w/ a column name 'stimulus'
#' @param timepoints Timepoints to use from the measured data

getMetaData = function(input_obj, stimulus, timepoints){
  require('Seurat')
  require('factoextra')
  require('matrixStats')
  require('cowplot')
  require('ggplot2')
  require('RColorBrewer')
  require('NbClust')
  require('devtools')
  require('NMF')
  require('gridExtra')
  require('Rfast')
  pca <- input_obj[['pca']]
  metadata <- as.data.frame(input_obj@meta.data)

  # Confine metadata to stimulus and Unstim 0 hr
  metadata <- metadata[metadata$stimulus == stimulus | metadata$timept == '0hr'| metadata$timept == '0.0hr',]

  # Add a numeric timept column to metadata
  metadata$timept_num = as.numeric(sapply(strsplit(metadata$timept,"h"), `[`, 1))
  metadata = metadata[trunc(metadata$timept_num, 5) %in% trunc(timepoints, 5),]

  return(metadata)
}
