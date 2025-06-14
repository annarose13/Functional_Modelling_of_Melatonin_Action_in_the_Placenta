### EXAMPLE HUMAN scRNAseq CORRESPONDENCE ANALYSIS ###

# This is an example script for the correspondence analysis conducted in the human data sets
# This was repeated for all Leiden Clusters in both the first and third trimester data.

library(pheatmap)
library(tibble)
library(dendextend)
library(gplots)
library(plyr)
library(dplyr)
library(BioQC)
library(forcats)
library(entropy)
library(ggplot2)


## Set the working directory:
setwd("~/RP2/Yang")

## Load in data
exp_data <- read.csv("yang_anndata.csvs/X.csv",header=F, sep = ',')
obs_data <-read.csv("yang_anndata.csvs/obs.csv",header=T)
var_data<- read.csv("yang_anndata.csvs/var.csv",header=T)

## Edit expression data by assigning row names as the cell ID column:
rownames(exp_data)<-obs_data$CellID
## Change cell type column from a character to a factor. This generates a new column in obs data called 'final_cell_type'
obs_data$final_cell_type<-as.factor(obs_data$unique.cell.type) 
## assign column names as genes names
colnames(exp_data)<-var_data$Gene
### DOUBLE CHECK ROW NAMES ARE SAMPLES/CELLS AND COLUMN NAMES ARE GENES BEFORE CONTINUING ###



####SUBSET DATA FOR EACH CLUSTER####
## Label the obs data by cell type/cluster:
END_1 <- obs_data[obs_data$final_cell_type == "END_1", ]
## Subset the expression data to include only expression data for the specific cluster (SHOULD CONTAIN ALL GENES AT THIS POINT)
END_1_X <- exp_data[na.omit(match(END_1$CellID,rownames(exp_data))),]
## sd scores for all clusters # IMPORTANT - Correlation will not work on genes with no variance - remove genes with SD of 0 (THIS REMOVES THE GENES)
sd.scores1<-apply(END_1_X,2,sd) #(2 = selects columns not rows)
END_1_X_sd<-END_1_X[,which(sd.scores1>0)]



#### CORRELATION MATRICIES ####
cor_data_END_1<-cor(END_1_X_sd,END_1_X_sd)
write.csv(cor_data_END_1, "~/RP2/Yang/Yang_Hypernetwork_Analysis/Correlation_matrices/cor_data_END_1.csv", row.names=FALSE)
#cor_data_END_1 <-read.csv("~/RP2/Yang/Yang_Hypernetwork_Analysis/Correlation_matrices/cor_data_END_1.csv")


####  LOAD IN MELATONIN GENES OF INTEREST ####
ALL_MELATONIN_Genes <- read.csv("~/RP2/melatonin_genes/MEL_Genes.csv", header=TRUE)
PORCINE_Genes_list <- ALL_MELATONIN_Genes$PORCINE_Genes
PORCINE_vars <- var_data[var_data$Gene %in% PORCINE_Genes_list, ]
## Subset the VT expression data to include only the expression data for melatonin related genes 
## This should be 42,215 rows (number of cells in Yang data) by the number of melatonin genes in selected list
PORCINE_X <- exp_data[, na.omit(match(PORCINE_vars$Gene,colnames(exp_data))),]

######END_1
## Subset the above expression data to only include cells labelled as END_1:
PORCINE_X_2 <- data.frame(CellID = rownames(PORCINE_X,), PORCINE_X) # Add a column of cell names to data frame
PORCINE_END_1_X <-subset(PORCINE_X_2, CellID %in%END_1$CellID) # Rows are allEND_1 cells, columns are still all melatonin genes of interest
PORCINE_END_1_X <- select(PORCINE_END_1_X, -CellID) # Remove cellID column again
# Calculate PORCINE_END_1 SD
sd.scores_PORCINE_X_END_1 <-apply(PORCINE_END_1_X,2,sd) # (2 = selects columns not rows)
PORCINE_X_END_1_sd<-PORCINE_END_1_X[,which(sd.scores_PORCINE_X_END_1>0)] # removes any with SD of 0




#### MELATONIN HYPERGRAPH ####

## Generates a matrix of Melatonin genes as rownames and allEND_1 genes (- the melatonin genes of interest) as column names
PORCINE_END_1_DEGs<- colnames(PORCINE_X_END_1_sd)
END_1_PORCINE_melatonin_result_df <- data.frame() #Create empty data frame
cor_data_END_1_PORCINE <- cor_data_END_1[, -na.omit(match(PORCINE_END_1_DEGs,colnames(cor_data_END_1)))] #This removes the melatonin genes (DEGs) from the colnames 
cor_data_END_1_PORCINE <- cor_data_END_1_PORCINE[na.omit(match(PORCINE_END_1_DEGs, rownames(cor_data_END_1))),] # This ensures only the melatonin genes of interest are in the rows
bin_END_1_PORCINE<-abs(cor_data_END_1_PORCINE) # Create correlation matrix w/out signs
bin_END_1_PORCINE[which(bin_END_1_PORCINE>sd(cor_data_END_1_PORCINE))]<-1 #Binarise the correlation matrix
bin_END_1_PORCINE[which(bin_END_1_PORCINE!=1)]<-0 # this creates the hypergraph incidence matrix
hyp_END_1_PORCINE <-bin_END_1_PORCINE %*% t(bin_END_1_PORCINE) #adjacency matrix 
END_1_PORCINE_melatonin_entropy_result<- entropy(hyp_END_1_PORCINE)
END_1_PORCINE_melatonin_result_df <- rbind(END_1_PORCINE_melatonin_result_df,END_1_PORCINE_melatonin_entropy_result)



#SAVE HYPERGRAPHS
#write.csv(hyp_END_1_PORCINE, "~/RP2/RP2/Yang/Yang_Hypernetwork_Analysis/Porcine_hypergraphs/hyp_END_1_PORCINE.csv", row.names=FALSE)
#hyp_END_1_PORCINE <- read.csv( "~/RP2/Yang/Yang_Hypernetwork_Analysis/Porcine_hypergraphs/hyp_END_1_PORCINE.csv")


#Plot adjacency matrix heatmap
#Use cutree_rows and cutree_columns for how many different clusters you want to split it into (2)
#DO THESE ONE AT A TIME TO CHECK THAT CORRECT CENTRAL CLUSTER HAS BEEN ANNOTATED!
### END_1
END_1_PORCINE_hm <- pheatmap::pheatmap(hyp_END_1_PORCINE, fontsize_row = 2,fontsize_col = 2,cutree_rows = 2,cutree_cols = 2, main ="Third Trimester END_1")
END_1_PORCINE_hm_clusters <- as.data.frame(cutree(END_1_PORCINE_hm$tree_row,2))
#optional creating column with the names of the nodes
END_1_PORCINE_hm_clusters$Gene <- rownames(END_1_PORCINE_hm_clusters)
#Separating out into 2 dataframes, each containing names of nodes in that cluster
#I just check the number of datapoints to see which cluster is which (normally pretty obvious from heatmap)
END_1_PORCINE_hm_cluster_1 <- as.data.frame(END_1_PORCINE_hm_clusters[which(END_1_PORCINE_hm_clusters$`cutree(END_1_PORCINE_hm$tree_row, 2)` == 2),])
END_1_PORCINE_hm_cluster_2 <- as.data.frame(END_1_PORCINE_hm_clusters[which(END_1_PORCINE_hm_clusters$`cutree(END_1_PORCINE_hm$tree_row, 2)` == 1),])


## Overall - this is a list of the genes in the central cluster
# Put this into a character vector
END_1_PORCINE_central_cluster<- END_1_PORCINE_hm_cluster_1$Gene
#Restrict the rows to only include genes in central cluster
row_restricted_END_1_PORCINE <- bin_END_1_PORCINE[na.omit(match(END_1_PORCINE_central_cluster,rownames(bin_END_1_PORCINE))),]
#Restrict columns to only include genes from END_1 that are correlated with 95% of the central genes
galois_95_END_1_PORCINE <- row_restricted_END_1_PORCINE[, which(colSums(row_restricted_END_1_PORCINE) >= 0.95 * length(END_1_PORCINE_central_cluster))]
galois_95_END_1_PORCINE <- as.data.frame(galois_95_END_1_PORCINE)
galois_95_END_1_PORCINE <- tibble::rownames_to_column(galois_95_END_1_PORCINE)

write.csv(galois_95_END_1_PORCINE, "~/RP2/Yang/Yang_Hypernetwork_Analysis/Yang_Galois_results/YANG_95_END_1_PORCINE_galois.csv", row.names=FALSE)

