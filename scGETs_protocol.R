# scGETs imputation method
# Ksheu, Nov 2024

#############################################################

# Quickstart----

#############################################################

library(devtools)
install_github("ksheu/scGETs")
library(scGETs)


setwd("F:///BACKUP_USB_20200710_active/Projects_writing/trajectory_method/")
macro = readRDS("./input/macrophage_M0_rep2only_500genes_DBEC_scGETs_LPS.rds")
select_timepoints = c(0.0, 0.25, 1, 3, 8)
metadata = getMetaData(macro, stimulus = "LPS", timepoints= select_timepoints)
undebug(scREALTIME)
reconst = scREALTIME(input_obj = macro, metadata = metadata, timepoints = select_timepoints, stimulus = "LPS",
                     num_archetypes = 250, num_sim_pts = 100, num_trajectories = 1000, 
                     reduction = 'pca', consensus_measure = 'median', interpolant = 'spline', 
                     data_norm = "ISnorm", prob_method = 'distance', distance_metric = 'euclidean' ,
                     varFilter = T, exp_prob = 1) 
reconstructed_pc = reconst$reconstructed_trajectories
#############################################################

# Before you begin ----

#############################################################


#############################################################
# 0.1. Set up enviroment -----
#############################################################

# Install and load packages
library(data.table);library(reshape2);library(pheatmap);library(resample) # data formatting packages
library(RColorBrewer);library(ggpubr);library(ggplot2);library(pheatmap) # visualization packages
library(Seurat) #single-cell analysis package


#############################################################
# 0.2. Load and subset the dataset of interest, format to Seurat Object -----
#############################################################

#pulling data from the prior Seurat object
if(0){
  setwd("F://scRNAseq_macro/scRNAseq_macro/")
  macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
  macro = subset(macro, subset= stimulus=="LPS"|timept=="0.0hr")
  gene = "Tnf"
  VlnPlot(object = subset(macro, subset= stimulus=="LPS"|timept=="0.0hr"), pt.size = 0.5, cols = rep("#00BA38",7),
          features = c(gene), group.by = "timept", assay = "ISnorm" ) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("LPS")
  
  my.counts.raw = macro@assays$RNA@counts 
  my.counts.ISnorm = macro@assays$ISnorm@counts 
  my.counts.ISnorm.data = macro@assays$ISnorm@data #log2+1
  metadata = macro@meta.data
  
  # saving to demonstrate formatting data into Seurat object
  setwd("F:///BACKUP_USB_20200710_active/Projects_writing/trajectory_method/")
  
  # write.table(as.matrix(my.counts.raw), "./input/my.counts_raw_LPS.txt", row.names = T,col.names = T,sep = "\t",quote = F )
  # write.table(as.matrix(my.counts.ISnorm), "./input/my.counts_ISnorm_LPS.txt", row.names = T,col.names = T,sep = "\t",quote = F )
  # write.table(as.matrix(my.counts.ISnorm.data), "./input/my.counts_ISnorm_LPS_log2.txt", row.names = T,col.names = T,sep = "\t",quote = F )
  # write.table(metadata, "./input/my.metadata_LPS.txt", row.names = T,col.names = T,sep = "\t",quote = F )
  
}


# Format data into Seurat object----

my.counts.raw = read.delim("./input/my.counts_raw_LPS.txt")
my.counts.ISnorm = read.delim("./input/my.counts_ISnorm_LPS.txt")
my.counts.ISnorm.data= read.delim("./input/my.counts_ISnorm_LPS_log2.txt")
metadata = read.delim("./input/my.metadata_LPS.txt")

head(my.counts.raw, 15)[1:5] 
head(my.counts.ISnorm, 15)[1:5] 
head(my.counts.ISnorm.data, 15)[1:5] 
head(metadata, 10) 

macro <- CreateSeuratObject(counts = my.counts.raw, project = "rhapsody")
macro[["ISnorm"]] <- CreateAssayObject(counts = my.counts.ISnorm)
macro[["ISnorm"]]@data = as.matrix(my.counts.ISnorm.data)
macro@meta.data = metadata

#run an unscaled, uncentered PCA
macro <- FindVariableFeatures(object = macro, assay = "ISnorm")
macro <- ScaleData(object = macro, assay = "ISnorm")
macro[["ISnorm"]]@scale.data = as.matrix(macro[["ISnorm"]]@data)
macro<- RunPCA(macro, assay = "ISnorm")

# saveRDS(macro,"./input/macrophage_M0_rep2only_500genes_DBEC_scGETs_LPS.rds")

gene = "Il1b"
VlnPlot(object = subset(macro), pt.size = 0.1, cols = rep("#00BA38",7),
        features = c(gene), group.by = "timept", assay = "ISnorm" ) +
  stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
  theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("LPS")
DimPlot(macro, reduction = "pca", group.by = "timept")
DimPlot(macro, reduction = "pca", group.by = "stimulus", split.by = "timept")

#############################################################
# 0.3. Pull metadata -----
#############################################################


# Get metadata of desired cells and subset data based on this----
macro = readRDS("./input/macrophage_M0_rep2only_500genes_DBEC_scGETs_LPS.rds")
select_timepoints = c(0.0, 0.25, 1, 3, 8)
timepoints= select_timepoints
stimulus = "LPS"

metadata <- as.data.frame(macro@meta.data)
metadata <- metadata[metadata$stimulus == stimulus | metadata$timept == '0hr'| metadata$timept == '0.0hr',]
metadata$timept_num = as.numeric(sapply(strsplit(metadata$timept,"h"), `[`, 1))

#############################################################
# 0.4. Retrieve matrix of cells of interest -----
#############################################################


# Retrieve the desired data subset (e.g. the data matrix with the timepoints to be included).
cells_by_timept <- list()
for(i in timepoints){
  index = paste("time_",i, "hr", sep = "")
  cells_by_timept[[index]] <- rownames(metadata)[metadata$timept_num == i]
} 

data = "ISnorm"
if(data == 'RNA'){
  RNA <- as.data.frame(macro@assays$RNA@data)
}else if(data == 'ISnorm'){
  RNA <- as.data.frame(macro@assays$ISnorm@data)
}

RNA <- t(RNA)
RNA <- RNA[rownames(RNA) %in% rownames(metadata),]


#############################################################

#Step-by-step Method Details----

#############################################################


#############################################################
# 1. Dimensionality Reduction via k-means and PCA -----
#############################################################

# Specify hyper parameters ----
input_obj = macro; metadata = metadata; timepoints = select_timepoints;
num_archetypes = 20; num_sim_pts = 100; num_trajectories = 1000;
data = "ISnorm"; reduction = 'pca'; stimulus = "LPS"; consensus_measure = 'median'; interpolant = 'spline';prob_method = 'distance'; distance_metric = 'euclidean' ; varFilter = T; exp_prob = 1

#k-means clustering ----
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
  
  if(i == timepoints[1]){
    cluster_counts <- as.data.frame(table(clusterings[[index]]$cluster))
  }else{
    cluster_counts <- cbind(cluster_counts, table(clusterings[[index]]$cluster))
  }
}

# Clean up cluster counts dataframe
col_omit <- c(1:length(timepoints))*2 - 1
cluster_counts <- cluster_counts[, -col_omit]
rownames(cluster_counts) <- c(1:num_archetypes)
rownames(cluster_counts) <- paste("bin", rownames(cluster_counts), sep="")
colnames(cluster_counts) <- timepoints
colnames(cluster_counts) <- paste("time_", colnames(cluster_counts), sep="")

# Calculate and store cluster densities from cluster_counts table
cluster_densities <- as.data.frame(prop.table(as.matrix(cluster_counts), 2))
cluster_densities <- cluster_densities^(exp_prob)

# Run an unscaled PCA and store within Seurat object
input_obj <- FindVariableFeatures(object = input_obj, assay = "ISnorm")
input_obj <- ScaleData(object = input_obj, assay = "ISnorm")
input_obj[["ISnorm"]]@scale.data = as.matrix(input_obj[["ISnorm"]]@data)
input_obj<- RunPCA(input_obj, assay = "ISnorm")

# Retrieve PCA (or other dimensionality reduction) results from Seurat object
if(toupper(reduction) == 'PCA'){
  pcscores = input_obj[['pca']]@cell.embeddings
}else if(toupper(reduction) == 'NMF'){
  pcscores = input_obj[['NMF']]@cell.embeddings
}else if(toupper(reduction) == 'ICA'){
  pcscores = input_obj[['ica']]@cell.embeddings
}

# DimPlot(input_obj, reduction = "pca", group.by = "timept")
# DimPlot(input_obj, reduction = "pca", group.by = "stimulus", split.by = "timept")

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
pcscores_stim$bin <- as.factor(cell_cluster_df$cell_cluster_df[match(rownames(pcscores_stim), rownames(cell_cluster_df))])
pcscores_stim$timebin_tag <- paste(as.numeric(pcscores_stim$timept), pcscores_stim$bin, sep = "_")

#plot intermediate
library(plyr)
ggplot(pcscores_stim, aes(PC_1, PC_2, color = as.factor(timept)))+geom_point(size=0.5,alpha=0.2)+
  theme_bw(base_size = 12)
# ggplot(pcscores_stim, aes(PC_1, PC_2, color = as.factor(bin)))+geom_point(size=0.5,alpha=0.5)+
#   theme_bw(base_size = 12)+facet_wrap(~timept, scales = "free")+
#   stat_ellipse(level = 0.9, geom="polygon",aes(fill=bin), alpha=0.1)  
#   #geom_polygon(data = hulls, alpha = 0.5) #+scale_color_distiller(palette = "Spectral")	

if(consensus_measure == 'mean'){
  aggreg_pcscores <- aggregate(pcscores_stim[, 1:(ncol(pcscores_stim)-3)], list(pcscores_stim$timebin_tag), mean)
}else if(consensus_measure == 'median'){
  aggreg_pcscores <- aggregate(pcscores_stim[, 1:(ncol(pcscores_stim)-3)], list(pcscores_stim$timebin_tag), median)
}
aggreg_pcscores$timept <- sapply(strsplit(aggreg_pcscores$Group.1, split = "_"), `[`, 1)
aggreg_pcscores$bin <- sapply(strsplit(aggreg_pcscores$Group.1, split = "_"), `[`, 2)
col_orders = c(1,(ncol(aggreg_pcscores)), (ncol(aggreg_pcscores)-1), 2:(ncol(aggreg_pcscores)-2))
aggreg_pcscores <- aggreg_pcscores[order(as.numeric(aggreg_pcscores$bin)),]
aggreg_pcscores <- aggreg_pcscores[order(aggreg_pcscores$timept),col_orders]
rownames(aggreg_pcscores) = aggreg_pcscores$Group.1


#plot intermediate
ggplot(aggreg_pcscores, aes(PC_1, PC_2, color = as.factor(timept)))+geom_point(alpha=0.5)+
  theme_bw(base_size = 12)
  
#############################################################
# 2. Transition Probability Matrices -----
#############################################################

#calulation of transition probability matrices-----
distances = as.matrix(dist(aggreg_pcscores[,4:ncol(aggreg_pcscores)],method = distance_metric))
rownames(distances) = aggreg_pcscores$Group.1
colnames(distances) = aggreg_pcscores$Group.1

library(pheatmap)
pheatmap(distances, cluster_rows = F, cluster_cols = F,
         annotation_row = aggreg_pcscores[,c(3), drop=F],
         annotation_col = aggreg_pcscores[,c(3), drop=F],
         show_rownames = F, show_colnames = F)

walks <- list()
walk_probs <- matrix(nrow = num_trajectories, ncol = 2)
walk_probs[,1] = c(1:num_trajectories)


if(prob_method == 'distance'){
  print("Generating random walks based on distance between clusters")
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
        
        # Normalize so larger distances have lower probability
        probs = (1/probs)^exp_prob
        next_value = sample(c(1:num_archetypes), size = 1, prob = probs)
      }
      
      path <- c(path, next_value)
    }
    walks[[i]] <- path
    
  }
}


# Checking number of unique walks
num_unique_traj = length(unique(( walks)))
print(paste(length(unique(( walks))), "unique walks/trajectories"))

# Make table of unique walks for us to keep track of
walk_frequencies <- as.data.frame(table(sapply( walks, paste0, collapse="")))
colnames(walk_frequencies)[1] = 'Walk'

#############################################################
# 3. Spline Fitting -----
#############################################################

#random walks specify which values from the aggreg_pcscores matrix to use as the anchor points for spline fitting

spline_pts <- list()
for(pc in 1:(ncol(aggreg_pcscores)-3)){
  id = paste("pc", pc, sep="_")
  spline_pts[[id]] <- matrix(nrow = num_trajectories, ncol = length(timepoints))
  #View(spline_pts)
  for(tr in 1:num_trajectories){
    path = walks[[tr]]
    for(j in 1:length(timepoints)){
      time = timepoints[j]
      bin = path[j]
      spline_pts[[id]][tr, j] <- aggreg_pcscores[round(as.numeric(aggreg_pcscores$timept),3) == round(time,3) & aggreg_pcscores$bin == bin, pc+3]
    }
  }
}

# Now, we get number of unique paths * number of trajectories * 50 PCs
# num_splines =dim(walk_frequencies)[1] * 50
# print(paste("Number of splines: ", dim(walk_frequencies)[1] * 50))
dup_rows = duplicated(spline_pts[[1]])
walk_probs = walk_probs[!dup_rows, ]

for(i in 1:length(spline_pts)){
  spline_pts[[i]] <- spline_pts[[i]][!dup_rows,]
  
}

# Add the timepts measured to num_sim_pts if not already included

sim_times = seq(min(timepoints),max(timepoints),length.out = num_sim_pts)
sim_times
for(i in timepoints){ 
  if(!(i %in% sim_times)){
    sim_times <- c(sim_times, i)
    num_sim_pts = num_sim_pts + 1
  }
}
sim_times <- sort(sim_times)

simulated = matrix(ncol = (ncol(aggreg_pcscores)-3), nrow = num_sim_pts*dim(walk_frequencies)[1])
simulation_mat_time = rep(sim_times, dim(walk_frequencies)[1])

#do spline interpolation
for(i in 1:(ncol(aggreg_pcscores)-3)){ # Loop over PCs
  preds = c()
  for(walk in 1:nrow(spline_pts[[i]])){
    if(interpolant == 'spline'){
      spline <- smooth.spline(timepoints, spline_pts[[i]][walk,], cv = T, all.knots = T) 
      preds <- c(preds, predict(spline, sim_times)$y)
    } else {
      print('Unspecified interpolant')
    }
  }
  simulated[,i] = preds
}

#plot fitted splines
if(0){
  data = input_obj@meta.data
  pca = data.frame(input_obj@reductions$pca@cell.embeddings)
  data= cbind(data, pca[ match(rownames(data), rownames(input_obj@reductions$pca)),])
  
  simulated.annot = data.frame(simulated)
  simulated.annot$time <- simulation_mat_time
  simulated.annot$path <- rep(1:dim(walk_frequencies)[1], each = num_sim_pts)
  
  ggplot(NULL)+
    geom_point(data = data[grepl("LPS|Unstim",data$stimulus),], mapping = aes(PC_1,PC_2, color = as.factor(timept)),size=1)+
    geom_path(data=simulated.annot, mapping=aes(X1,X2,group=path),alpha=0.05)+
    theme_bw(base_size = 12)+xlab("PC1")+ylab("PC2")

  
  ggplot(NULL)+
    geom_point(data = aggreg_pcscores, mapping = aes(PC_1,PC_2, color = as.factor(timept)),size=1)+
    geom_path(data=simulated.annot, mapping=aes(X1,X2,group=path),alpha=0.05)+
    theme_bw(base_size = 12)+xlab("PC1")+ylab("PC2")
}

#############################################################
# 4. Recovering gene expression trajectories -----
#############################################################

# Deconvolute to gene expression
loadings <- input_obj[[reduction]]@feature.loadings

reconstructed_pc <- t(loadings %*% t(simulated))
reconstructed_pc <- as.data.frame(reconstructed_pc)
reconstructed_pc$time <- simulation_mat_time
reconstructed_pc$path <- rep(1:dim(walk_frequencies)[1], each = num_sim_pts)

#plot output---------------------

reconstructed_pc.traj = (reconstructed_pc)
reconstructed_pc.traj$stimulus = stimulus
reconstructed_pc.traj$type = "M0"
reconstructed_pc.traj$path_stim = paste0(reconstructed_pc.traj$path,"_", reconstructed_pc.traj$stimulus)

mat.numbers = reconstructed_pc.traj[,!grepl("time|path|stimulus|type|path_stim", colnames(reconstructed_pc.traj))]
mat.numbers = apply(mat.numbers, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))) #rescale each gene column 0-1 over all stims
#mat.numbers[mat.numbers<=0] <-0 #set negative to 0
mat.numbers = cbind(mat.numbers, reconstructed_pc.traj[,grepl("time|path|stimulus|type|path_stim", colnames(reconstructed_pc.traj))])

gene = "Cxcl10"
ggplot(mat.numbers[grepl("",mat.numbers$stimulus)&grepl("",mat.numbers$type),], 
       aes(time,get(gene), group = path_stim)) +
  geom_vline(xintercept = c(0,0.25,1,3,8), linetype="dotted")+ #theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  # geom_vline(xintercept = c(0,0.5,1,3,5,8), linetype="dotted")+
  geom_line(aes(group = as.factor(path_stim), color = type), alpha = 0.05)+
  theme_bw(base_size = 12)+theme(legend.position = "none")+ylab(gene)



#plot trajectories in heatmap form
gene = "Ccl5"
mat.numbers.dcast = dcast(mat.numbers, stimulus+path~time, value.var = gene)
rownames(mat.numbers.dcast) = paste0(mat.numbers.dcast$path,"_", mat.numbers.dcast$stimulus)
colors_list = list(stimulus = c(CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"))
annot.frame = mat.numbers.dcast[,c(1,2)] #[,c(4,5,6,8,13),drop=F] #c(13,4,5,6,8,9,10)
mat.numbers.dcast = mat.numbers.dcast[grepl("", mat.numbers.dcast$stimulus),]

pheatmap(as.matrix(mat.numbers.dcast[,-c(1,2)]),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         breaks = c(0,seq(0.01,0.99,length=100),1),
         cluster_cols = F, cluster_rows = T,
         show_colnames = F, show_rownames = F,
         clustering_method = "ward.D2",
         main = gene, annotation_row = annot.frame,
         annotation_colors = colors_list)

#plot trajectories heatmap, multiple genes
gene = "Cxcl10"
X2.mat.numbers.dcast = dcast(mat.numbers, stimulus+path~time, value.var = gene)
rownames(X2.mat.numbers.dcast) = paste0(X2.mat.numbers.dcast$path,"_", X2.mat.numbers.dcast$stimulus)

gene = "Il1b"#"Gna15","Egr3"
X3.mat.numbers.dcast = dcast(mat.numbers, stimulus+path~time, value.var = gene)
rownames(X3.mat.numbers.dcast) = paste0(X3.mat.numbers.dcast$path,"_", X3.mat.numbers.dcast$stimulus)

gene = "Ccl5"#"Gna15","Egr3"
X4.mat.numbers.dcast = dcast(mat.numbers, stimulus+path~time, value.var = gene)
rownames(X4.mat.numbers.dcast) = paste0(X4.mat.numbers.dcast$path,"_", X4.mat.numbers.dcast$stimulus)

gene = "Nfkbia"
mat.numbers.dcast = dcast(mat.numbers, stimulus+path~time, value.var = gene)
rownames(mat.numbers.dcast) = paste0(mat.numbers.dcast$path,"_", mat.numbers.dcast$stimulus)

mat.numbers2.dcast = cbind(X2.mat.numbers.dcast,
                           X3.mat.numbers.dcast[,-c(1,2)],
                           X4.mat.numbers.dcast[,-c(1,2)],
                           mat.numbers.dcast[,-c(1,2)])
colors_list = list(path_stim = c(CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"))
# colors_list = list(path_stim = c(`0_rep2only_TNF`="darkred", `1_IFNg_TNF`="#00BA38", `2_IL4_gt80_TNF`="#619CFF"))
annot.frame = mat.numbers.dcast[,c(2),drop=F] #c(13,4,5,6,8,9,10)
annot.frame$path_stim = gsub("_..*","", annot.frame$path_stim)

pheatmap(as.matrix(mat.numbers2.dcast[,-c(1,2)]), scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         breaks = c(0,seq(0.01,0.99,length=100),1),
         cluster_cols = F, cluster_rows = T,
         show_colnames = F, show_rownames = F,
         clustering_method = "ward.D2",
         gaps_col = c(103,206,309),
         main = gene, annotation_row = annot.frame,
         annotation_colors = colors_list)

###################################################
# Optional: Calculate trajectory features ----
###################################################
#Calculate peak induction of each gene
mat.numbers = reconstructed_pc[,!grepl("time|path", colnames(reconstructed_pc))]
mat.numbers = apply(mat.numbers, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))) #rescale each gene column 0-1 over all stims
mat.numbers = cbind(mat.numbers, reconstructed_pc[,grepl("time|path", colnames(reconstructed_pc))])

dynamics = data.frame()
for (i in colnames(mat.numbers)[!grepl("time|path", colnames(mat.numbers))]){
  print(i)
  gene_name = i
  A.subset = reshape2::dcast(mat.numbers, path~time, value.var = i)
  my.dataframe = cbind(label = stimulus, A.subset[,-1])
  peak_amp <- apply(my.dataframe[,-1], 1, max) 
  tmp = data.frame(peak_amp =peak_amp, stimulus =stimulus, gene = gene_name)
  dynamics <- rbind(dynamics, tmp)
}

gene = "Cxcl10"
ggplot(dynamics[grepl(paste0(gene,"$"),dynamics$gene),], aes(stimulus, peak_amp))+
  geom_violin(aes(color = stimulus))+geom_point(position = "jitter",alpha = 0.5)+
  stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
  theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")

#peak fold change
dynamics = data.frame()
for (i in colnames(mat.numbers)[!grepl("time|path", colnames(mat.numbers))]){
  print(i)
  gene_name = i
  A.subset = reshape2::dcast(mat.numbers, path~time, value.var = i)
  my.dataframe = cbind(label = stimulus, A.subset[,-1])
  peak_amp <- apply(my.dataframe[,-1], 1, max) 
  tmp = data.frame(peak_amp =peak_amp, 
                   peak_amp_lfc = log2((peak_amp/(my.dataframe[,2]+0.01))+1), #2nd col = time0
                   peak_amp_fc = (peak_amp/(my.dataframe[,2]+0.01)),
                   time0_amp = my.dataframe[,2],
                   stimulus =stimulus, gene = gene_name) 
  dynamics <- rbind(dynamics, tmp)
}

gene = "Tnf"
ggplot(dynamics[grepl(paste0(gene,"$"),dynamics$gene),], aes(stimulus, peak_amp_lfc))+
  geom_violin(aes(color = stimulus))+geom_point(position = "jitter", alpha = 0.5)+
  stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
  theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")


#Speed at time 1hr
dynamics = data.frame()
for (i in colnames(mat.numbers)[!grepl("time|path", colnames(mat.numbers))]){
  print(i)
  gene_name = i
  A.subset = reshape2::dcast(mat.numbers, path~time, value.var = i)
  my.dataframe = cbind(label = stimulus, A.subset[,-1])
  timept_tangent = which(colnames(my.dataframe[,-1]) == "1") #for slope at 1hr

  timeseg <- as.numeric(names(my.dataframe[,-1])[round(timept_tangent)+2]) - 
    as.numeric(names(my.dataframe[,-1])[round(timept_tangent)-2])
  rise <- my.dataframe[,-1][,round(timept_tangent)+2]- my.dataframe[,round(timept_tangent)-2]
  
  tmp = data.frame(speed1hr = (rise/timeseg), stimulus =stimulus, gene = gene_name) 
  dynamics <- rbind(dynamics, tmp)
}

gene = "Tnf"
ggplot(dynamics[grepl(paste0(gene,"$"),dynamics$gene),], aes(stimulus, speed1hr))+
  geom_violin(aes(color = stimulus))+geom_point(position = "jitter",alpha=0.4)+
  stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
  theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")


# Integral, total mRNA
dynamics = data.frame()
for (i in colnames(mat.numbers)[!grepl("time|path", colnames(mat.numbers))]){
  tryCatch(
    expr = {
      print(i)
      gene_name = i
      A.subset = reshape2::dcast(mat.numbers, path~time, value.var = i)
      my.dataframe = cbind(label = stimulus, A.subset[,-1])
      
      time <- as.numeric(colnames(my.dataframe[-1]))
      integral <- apply(my.dataframe[,-1], 1, function(x) unlist(integrate(approxfun(time, x), range(time)[1], range(time)[2],rel.tol =.Machine$double.eps^.2))$value)
      
      tmp = data.frame(integral = integral, stimulus =stimulus, gene = gene_name) 
      dynamics <- rbind(dynamics, tmp)
    }, error = function(e){
      message('Caught an error!')
      integral <- NA
      tmp = data.frame(integral = integral, stimulus =stimulus, gene = gene_name) 
      dynamics <- rbind(dynamics, tmp)
    }
  )
}
gene = "Tnf"
ggplot(dynamics[grepl(paste0(gene,"$"),dynamics$gene),], aes(stimulus, integral))+
  geom_violin(aes(color = stimulus))+geom_point(position = "jitter", alpha = 0.5)+
  stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
  theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")



