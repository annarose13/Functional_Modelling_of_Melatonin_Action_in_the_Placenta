#### Example Hypergraph Analysis of Human scRNA-seq data ####

# This example shows workflow for the CTB and STB clusters in the third trimester data set.
# The process for all other cell types in this data set and all data from the first trimester was the same.

library(dendextend)
library(gplots)
library(plyr)
library(dplyr)
library(BioQC)
library(forcats)
library(entropy)
library(ggplot2)


## Set the working directory:
setwd()

## Load in data (these files are generated in First_Trimester_scRNAseq_Analysis.ipynb)
# The third trimester data files can be generated in Third_Trimester_scRNAseq_Analysis.ipynb 
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



        ###SUBSET DATA FOR EACH CLUSTER
## Label the obs data by cell type/cluster:
CTB_1 <- obs_data[obs_data$final_cell_type == "CTB_1", ]
## Subset the expression data to include only expression data for the specific cluster (SHOULD CONTAIN ALL GENES AT THIS POINT)
CTB_1_X <- exp_data[na.omit(match(CTB_1$CellID,rownames(exp_data))),]
## sd scores for all clusters # IMPORTANT - Correlation will not work on genes with no variance - remove genes with SD of 0 (THIS REMOVES THE GENES)
sd.scores1<-apply(CTB_1_X,2,sd) #(2 = selects columns not rows)
CTB_1_X_sd<-CTB_1_X[,which(sd.scores1>0)]

## Label the obs data by cell type/cluster:
CTB_2 <- obs_data[obs_data$final_cell_type == "CTB_2", ]
## Subset the expression data to include only expression data for the specific cluster (SHOULD CONTAIN ALL GENES AT THIS POINT)
CTB_2_X <- exp_data[na.omit(match(CTB_2$CellID,rownames(exp_data))),]
## sd scores for all clusters # IMPORTANT - Correlation will not work on genes with no variance - remove genes with SD of 0 (THIS REMOVES THE GENES)
sd.scores1<-apply(CTB_2_X,2,sd) #(2 = selects columns not rows)
CTB_2_X_sd<-CTB_2_X[,which(sd.scores1>0)]

## Label the obs data by cell type/cluster:
CTB_3 <- obs_data[obs_data$final_cell_type == "CTB_3", ]
## Subset the expression data to include only expression data for the specific cluster (SHOULD CONTAIN ALL GENES AT THIS POINT)
CTB_3_X <- exp_data[na.omit(match(CTB_3$CellID,rownames(exp_data))),]
## sd scores for all clusters # IMPORTANT - Correlation will not work on genes with no variance - remove genes with SD of 0 (THIS REMOVES THE GENES)
sd.scores1<-apply(CTB_3_X,2,sd) #(2 = selects columns not rows)
CTB_3_X_sd<-CTB_3_X[,which(sd.scores1>0)]

## Label the obs data by cell type/cluster:
CTB_4 <- obs_data[obs_data$final_cell_type == "CTB_4", ]
## Subset the expression data to include only expression data for the specific cluster (SHOULD CONTAIN ALL GENES AT THIS POINT)
CTB_4_X <- exp_data[na.omit(match(CTB_4$CellID,rownames(exp_data))),]
## sd scores for all clusters # IMPORTANT - Correlation will not work on genes with no variance - remove genes with SD of 0 (THIS REMOVES THE GENES)
sd.scores1<-apply(CTB_4_X,2,sd) #(2 = selects columns not rows)
CTB_4_X_sd<-CTB_4_X[,which(sd.scores1>0)]

## Label the obs data by cell type/cluster:
STB_1 <- obs_data[obs_data$final_cell_type == "STB_1", ]
## Subset the expression data to include only expression data for the specific cluster (SHOULD CONTAIN ALL GENES AT THIS POINT)
STB_1_X <- exp_data[na.omit(match(STB_1$CellID,rownames(exp_data))),]
## sd scores for all clusters # IMPORTANT - Correlation will not work on genes with no variance - remove genes with SD of 0 (THIS REMOVES THE GENES)
sd.scores1<-apply(STB_1_X,2,sd) #(2 = selects columns not rows)
STB_1_X_sd<-STB_1_X[,which(sd.scores1>0)]


        ### CORRELATION MATRICIES ###
cor_data_CTB_1<-cor(CTB_1_X_sd,CTB_1_X_sd) # all genes in cluster, this should be 4568x4568 genes for first cluster
cor_data_CTB_2<-cor(CTB_2_X_sd,CTB_2_X_sd)
cor_data_CTB_3<-cor(CTB_3_X_sd,CTB_3_X_sd)
cor_data_CTB_4<-cor(CTB_4_X_sd,CTB_4_X_sd)
cor_data_STB_1<-cor(STB_1_X_sd,STB_1_X_sd)

#cor_data_CTB_1 <- read.csv("Yang_Hypernetwork_Analysis/Correlation_matrices/cor_data_CTB_1.csv")
#cor_data_CTB_2 <- read.csv("Yang_Hypernetwork_Analysis/Correlation_matrices/cor_data_CTB_2.csv")
#cor_data_CTB_3 <- read.csv("Yang_Hypernetwork_Analysis/Correlation_matrices/cor_data_CTB_3.csv")
#cor_data_CTB_4 <- read.csv("Yang_Hypernetwork_Analysis/Correlation_matrices/cor_data_CTB_4.csv")
#cor_data_STB_1 <- read.csv("Yang_Hypernetwork_Analysis/Correlation_matrices/cor_data_STB_1.csv")


        ###  LOAD IN MELATONIN GENES OF INTEREST ###
ALL_MELATONIN_Genes <- read.csv(paste0(getwd(),"List_of_melatonin_related_genes.csv", header=TRUE) #This file is in the repo

CONTROL_Genes_list <- ALL_MELATONIN_Genes$CONTROL_Genes
CONTROL_vars <- var_data[var_data$Gene %in% CONTROL_Genes_list, ]
## Subset the VT expression data to include only the expression data for melatonin related genes 
## This should be 42,215 rows (number of cells in Yang data) by the number of melatonin genes in selected list
CONTROL_X <- exp_data[, na.omit(match(CONTROL_vars$Gene,colnames(exp_data))),]
    
######CTB_1
## Subset the above expression data to only include cells labelled as CTB_1:
CONTROL_X_2 <- data.frame(CellID = rownames(CONTROL_X,), CONTROL_X) # Add a column of cell names to data frame
CONTROL_CTB_1_X <-subset(CONTROL_X_2, CellID %in%CTB_1$CellID) # Rows are allCTB_1 cells, columns are still all melatonin genes of interest
CONTROL_CTB_1_X <- select(CONTROL_CTB_1_X, -CellID) # Remove cellID column again
# Calculate CONTROL_CTB_1 SD
sd.scores_CONTROL_X_CTB_1 <-apply(CONTROL_CTB_1_X,2,sd) # (2 = selects columns not rows)
CONTROL_X_CTB_1_sd<-CONTROL_CTB_1_X[,which(sd.scores_CONTROL_X_CTB_1>0)] # removes any with SD of 0

######CTB_2
## Subset the above expression data to only include cells labelled as CTB_2:
CONTROL_X_2 <- data.frame(CellID = rownames(CONTROL_X,), CONTROL_X) # Add a column of cell names to data frame
CONTROL_CTB_2_X <-subset(CONTROL_X_2, CellID %in%CTB_2$CellID) # Rows are allCTB_2 cells, columns are still all melatonin genes of interest
CONTROL_CTB_2_X <- select(CONTROL_CTB_2_X, -CellID) # Remove cellID column again
# Calculate CONTROL_CTB_2 SD
sd.scores_CONTROL_X_CTB_2 <-apply(CONTROL_CTB_2_X,2,sd) # (2 = selects columns not rows)
CONTROL_X_CTB_2_sd<-CONTROL_CTB_2_X[,which(sd.scores_CONTROL_X_CTB_2>0)] # removes any with SD of 0

######CTB_3
## Subset the above expression data to only include cells labelled as CTB_3:
CONTROL_X_2 <- data.frame(CellID = rownames(CONTROL_X,), CONTROL_X) # Add a column of cell names to data frame
CONTROL_CTB_3_X <-subset(CONTROL_X_2, CellID %in%CTB_3$CellID) # Rows are allCTB_3 cells, columns are still all melatonin genes of interest
CONTROL_CTB_3_X <- select(CONTROL_CTB_3_X, -CellID) # Remove cellID column again
# Calculate CONTROL_CTB_3 SD
sd.scores_CONTROL_X_CTB_3 <-apply(CONTROL_CTB_3_X,2,sd) # (2 = selects columns not rows)
CONTROL_X_CTB_3_sd<-CONTROL_CTB_3_X[,which(sd.scores_CONTROL_X_CTB_3>0)] # removes any with SD of 0

######CTB_4
## Subset the above expression data to only include cells labelled as CTB_4:
CONTROL_X_2 <- data.frame(CellID = rownames(CONTROL_X,), CONTROL_X) # Add a column of cell names to data frame
CONTROL_CTB_4_X <-subset(CONTROL_X_2, CellID %in%CTB_4$CellID) # Rows are allCTB_4 cells, columns are still all melatonin genes of interest
CONTROL_CTB_4_X <- select(CONTROL_CTB_4_X, -CellID) # Remove cellID column again
# Calculate CONTROL_CTB_4 SD
sd.scores_CONTROL_X_CTB_4 <-apply(CONTROL_CTB_4_X,2,sd) # (2 = selects columns not rows)
CONTROL_X_CTB_4_sd<-CONTROL_CTB_4_X[,which(sd.scores_CONTROL_X_CTB_4>0)] # removes any with SD of 0

######STB_1
## Subset the above expression data to only include cells labelled as STB_1:
CONTROL_X_2 <- data.frame(CellID = rownames(CONTROL_X,), CONTROL_X) # Add a column of cell names to data frame
CONTROL_STB_1_X <-subset(CONTROL_X_2, CellID %in%STB_1$CellID) # Rows are allSTB_1 cells, columns are still all melatonin genes of interest
CONTROL_STB_1_X <- select(CONTROL_STB_1_X, -CellID) # Remove cellID column again
# Calculate CONTROL_STB_1 SD
sd.scores_CONTROL_X_STB_1 <-apply(CONTROL_STB_1_X,2,sd) # (2 = selects columns not rows)
CONTROL_X_STB_1_sd<-CONTROL_STB_1_X[,which(sd.scores_CONTROL_X_STB_1>0)] # removes any with SD of 0






        ### HYPERGRAPHS ###

### CHECK EVERY SAMPLE AND SIZE!!!! ####
### sample = cluster_X_sd vars = total number of genes in cluster, size = mellist_x_cluster_sd vars = no. of mel genes in that cluster


###CTB_1
CTB_1_CONTROL_result_df <- data.frame() #Create empty data frame

for (i in 1:1000) {  # Outer loop (runs 1000 times)
  print(i)
  CTB_1_DEGs<- colnames(CTB_1_X_sd)[sample(1:4568,size=117)] # Randomly select genes from rows and label them DEGs
  cor_data_CTB_1_2 <- cor_data_CTB_1[, -na.omit(match(CTB_1_DEGs,colnames(cor_data_CTB_1)))] # This removes the random genes (DEGs) chosen in line above from the colnames
  cor_data_CTB_1_2 <- cor_data_CTB_1_2[na.omit(match(CTB_1_DEGs, rownames(cor_data_CTB_1))),] # This only keeps the random DEGs in the rownames (overall n obs/entries by 3577 columns)
  bin_CTB_1<-abs(cor_data_CTB_1_2) # Create correlation matrix w/out signs
  bin_CTB_1[which(bin_CTB_1>sd(cor_data_CTB_1_2))]<-1 #Binarise the correlation matrix
  bin_CTB_1[which(bin_CTB_1!=1)]<-0 # this creates the hypergraph incidence matrix
  hyp_CTB_1<-bin_CTB_1 %*% t(bin_CTB_1) #adjacency matrix
  CTB_1_entropy_result<- entropy(hyp_CTB_1)
  CTB_1_CONTROL_result_df <- rbind(CTB_1_CONTROL_result_df,CTB_1_entropy_result)
}
write.csv(CTB_1_CONTROL_result_df, paste0(getwd(), "/Yang_Hypernetwork_Analysis/CTB_1_CONTROL_entropy_results.csv"), row.names=FALSE)

      ### MELATONIN HYPERGRAPH ###

## Generates a matrix of Melatonin genes as rownames and allCTB_1 genes (- the melatonin genes of interest) as column names
CONTROL_CTB_1_DEGs<- colnames(CONTROL_X_CTB_1_sd)
CTB_1_CONTROL_melatonin_result_df <- data.frame() #Create empty data frame
cor_data_CTB_1_CONTROL <- cor_data_CTB_1[, -na.omit(match(CONTROL_CTB_1_DEGs,colnames(cor_data_CTB_1)))] #This removes the melatonin genes (DEGs) from the colnames 
cor_data_CTB_1_CONTROL <- cor_data_CTB_1_CONTROL[na.omit(match(CONTROL_CTB_1_DEGs, rownames(cor_data_CTB_1))),] # This ensures only the melatonin genes of interest are in the rows
bin_CTB_1_CONTROL<-abs(cor_data_CTB_1_CONTROL) # Create correlation matrix w/out signs
bin_CTB_1_CONTROL[which(bin_CTB_1_CONTROL>sd(cor_data_CTB_1_CONTROL))]<-1 #Binarise the correlation matrix
bin_CTB_1_CONTROL[which(bin_CTB_1_CONTROL!=1)]<-0 # this creates the hypergraph incidence matrix
hyp_CTB_1_CONTROL <-bin_CTB_1_CONTROL %*% t(bin_CTB_1_CONTROL) #adjacency matrix 
CTB_1_CONTROL_melatonin_entropy_result<- entropy(hyp_CTB_1_CONTROL)
CTB_1_CONTROL_melatonin_result_df <- rbind(CTB_1_CONTROL_melatonin_result_df,CTB_1_CONTROL_melatonin_entropy_result)

##PLOT MELATONIN HYPERGRAPH
#library(gplots)
#hm_CTB_1_CONTROL<-heatmap.2(hyp_CTB_1_CONTROL,trace="none",labRow = F,labCol=F)
#hm_CTB_1_CONTROL <- heatmap.2(hyp_CTB_1_CONTROL, trace = "none", labRow = FALSE, labCol = FALSE)
#dend_CTB_1_CONTROL<-as.hclust(hm_CTB_1_CONTROL$rowDendrogram)

### Z-scores ###
#Calculate z-scores for first hypergraphs
colnames(CTB_1_CONTROL_result_df) <- "CTB_1_CONTROL_Raw_Entropy_Values"
CTB_1_CONTROL_entropy_mean <- mean(CTB_1_CONTROL_result_df$CTB_1_CONTROL_Raw_Entropy_Values)
CTB_1_CONTROL_entropy_sd <- sd(CTB_1_CONTROL_result_df$CTB_1_CONTROL_Raw_Entropy_Values)
CTB_1_CONTROL_entropy_Z_scores <- (CTB_1_CONTROL_result_df$CTB_1_CONTROL_Raw_Entropy_Values -CTB_1_CONTROL_entropy_mean)/CTB_1_CONTROL_entropy_sd
# Add to column in entropy data frame:
z_scores_df <- data.frame(z_scores =CTB_1_CONTROL_entropy_Z_scores)
CTB_1_CONTROL_entropy_results <- cbind(CTB_1_CONTROL_result_df, z_scores_df)

### Calculate Z-score from the melatonin gene hypergraph:
CTB_1_CONTROL_MELATONIN_entropy_Z_score <- (CTB_1_CONTROL_melatonin_entropy_result-CTB_1_CONTROL_entropy_mean)/CTB_1_CONTROL_entropy_sd
### Calculate p-value and add to results file:
CTB_1_CONTROL_entropy_p_value <- 2*pnorm(-abs(CTB_1_CONTROL_MELATONIN_entropy_Z_score)) # 2* as it is two sided
#Add these to the results data frame:
CTB_1_CONTROL_entropy_results$CTB_1_CONTROL_melatonin_entropy_result <- NA
CTB_1_CONTROL_entropy_results$CTB_1_CONTROL_melatonin_entropy_result[1] <-CTB_1_CONTROL_melatonin_entropy_result
CTB_1_CONTROL_entropy_results$CONTROL_MELATONIN_entropy_Z_score <- NA
CTB_1_CONTROL_entropy_results$CONTROL_MELATONIN_entropy_Z_score[1] <-CTB_1_CONTROL_MELATONIN_entropy_Z_score
CTB_1_CONTROL_entropy_results$CTB_1_CONTROL_entropy_p_value <- NA
CTB_1_CONTROL_entropy_results$CTB_1_CONTROL_entropy_p_value[1] <-CTB_1_CONTROL_entropy_p_value
# Overwrite results file:
write.csv(CTB_1_CONTROL_entropy_results, paste0(getwd(), "/Yang_Hypernetwork_Analysis/CTB_1_CONTROL_entropy_results.csv"), row.names=FALSE)

### Plot the distribution:
hist(CTB_1_CONTROL_entropy_results$z_scores, breaks = 20, freq = FALSE, cex.lab=2, cex.axis=1.4, col = "skyblue", main = "Distribution of Z-Scores - CTB_1 CONTROL Third Trimester", xlab = "Z-Scores")
lines(density(CTB_1_CONTROL_entropy_results$z_scores), col = "red", lwd = 2)
#Add a line showing Z-score of melatonin:
abline(v =CTB_1_CONTROL_MELATONIN_entropy_Z_score, col = "green", lwd = 2)







###CTB_2
CTB_2_CONTROL_result_df <- data.frame() #Create empty data frame

for (i in 1:1000) {  # Outer loop (runs 1000 times)
  print(i)
  CTB_2_DEGs<- colnames(CTB_2_X_sd)[sample(1:4429,size=111)] # Randomly select genes from rows and label them DEGs
  cor_data_CTB_2_2 <- cor_data_CTB_2[, -na.omit(match(CTB_2_DEGs,colnames(cor_data_CTB_2)))] # This removes the random genes (DEGs) chosen in line above from the colnames
  cor_data_CTB_2_2 <- cor_data_CTB_2_2[na.omit(match(CTB_2_DEGs, rownames(cor_data_CTB_2))),] # This only keeps the random DEGs in the rownames (overall n obs/entries by 3577 columns)
  bin_CTB_2<-abs(cor_data_CTB_2_2) # Create correlation matrix w/out signs
  bin_CTB_2[which(bin_CTB_2>sd(cor_data_CTB_2_2))]<-1 #Binarise the correlation matrix
  bin_CTB_2[which(bin_CTB_2!=1)]<-0 # this creates the hypergraph incidence matrix
  hyp_CTB_2<-bin_CTB_2 %*% t(bin_CTB_2) #adjacency matrix
  CTB_2_entropy_result<- entropy(hyp_CTB_2)
  CTB_2_CONTROL_result_df <- rbind(CTB_2_CONTROL_result_df,CTB_2_entropy_result)
}
write.csv(CTB_2_CONTROL_result_df, paste0(getwd(), "/Yang_Hypernetwork_Analysis/CTB_2_CONTROL_entropy_results.csv"), row.names=FALSE)

### MELATONIN HYPERGRAPH ###

## Generates a matrix of Melatonin genes as rownames and allCTB_2 genes (- the melatonin genes of interest) as column names
CONTROL_CTB_2_DEGs<- colnames(CONTROL_X_CTB_2_sd)
CTB_2_CONTROL_melatonin_result_df <- data.frame() #Create empty data frame
cor_data_CTB_2_CONTROL <- cor_data_CTB_2[, -na.omit(match(CONTROL_CTB_2_DEGs,colnames(cor_data_CTB_2)))] #This removes the melatonin genes (DEGs) from the colnames 
cor_data_CTB_2_CONTROL <- cor_data_CTB_2_CONTROL[na.omit(match(CONTROL_CTB_2_DEGs, rownames(cor_data_CTB_2))),] # This ensures only the melatonin genes of interest are in the rows
bin_CTB_2_CONTROL<-abs(cor_data_CTB_2_CONTROL) # Create correlation matrix w/out signs
bin_CTB_2_CONTROL[which(bin_CTB_2_CONTROL>sd(cor_data_CTB_2_CONTROL))]<-1 #Binarise the correlation matrix
bin_CTB_2_CONTROL[which(bin_CTB_2_CONTROL!=1)]<-0 # this creates the hypergraph incidence matrix
hyp_CTB_2_CONTROL <-bin_CTB_2_CONTROL %*% t(bin_CTB_2_CONTROL) #adjacency matrix 
CTB_2_CONTROL_melatonin_entropy_result<- entropy(hyp_CTB_2_CONTROL)
CTB_2_CONTROL_melatonin_result_df <- rbind(CTB_2_CONTROL_melatonin_result_df,CTB_2_CONTROL_melatonin_entropy_result)

##PLOT MELATONIN HYPERGRAPH
#library(gplots)
#hm_CTB_2_CONTROL<-heatmap.2(hyp_CTB_2_CONTROL,trace="none",labRow = F,labCol=F)
#hm_CTB_2_CONTROL <- heatmap.2(hyp_CTB_2_CONTROL, trace = "none", labRow = FALSE, labCol = FALSE)
#dend_CTB_2_CONTROL<-as.hclust(hm_CTB_2_CONTROL$rowDendrogram)

### Z-scores ###
#Calculate z-scores for first hypergraphs
colnames(CTB_2_CONTROL_result_df) <- "CTB_2_CONTROL_Raw_Entropy_Values"
CTB_2_CONTROL_entropy_mean <- mean(CTB_2_CONTROL_result_df$CTB_2_CONTROL_Raw_Entropy_Values)
CTB_2_CONTROL_entropy_sd <- sd(CTB_2_CONTROL_result_df$CTB_2_CONTROL_Raw_Entropy_Values)
CTB_2_CONTROL_entropy_Z_scores <- (CTB_2_CONTROL_result_df$CTB_2_CONTROL_Raw_Entropy_Values -CTB_2_CONTROL_entropy_mean)/CTB_2_CONTROL_entropy_sd
# Add to column in entropy data frame:
z_scores_df <- data.frame(z_scores =CTB_2_CONTROL_entropy_Z_scores)
CTB_2_CONTROL_entropy_results <- cbind(CTB_2_CONTROL_result_df, z_scores_df)

### Calculate Z-score from the melatonin gene hypergraph:
CTB_2_CONTROL_MELATONIN_entropy_Z_score <- (CTB_2_CONTROL_melatonin_entropy_result-CTB_2_CONTROL_entropy_mean)/CTB_2_CONTROL_entropy_sd
### Calculate p-value and add to results file:
CTB_2_CONTROL_entropy_p_value <- 2*pnorm(-abs(CTB_2_CONTROL_MELATONIN_entropy_Z_score)) # 2* as it is two sided
#Add these to the results data frame:
CTB_2_CONTROL_entropy_results$CTB_2_CONTROL_melatonin_entropy_result <- NA
CTB_2_CONTROL_entropy_results$CTB_2_CONTROL_melatonin_entropy_result[1] <-CTB_2_CONTROL_melatonin_entropy_result
CTB_2_CONTROL_entropy_results$CONTROL_MELATONIN_entropy_Z_score <- NA
CTB_2_CONTROL_entropy_results$CONTROL_MELATONIN_entropy_Z_score[1] <-CTB_2_CONTROL_MELATONIN_entropy_Z_score
CTB_2_CONTROL_entropy_results$CTB_2_CONTROL_entropy_p_value <- NA
CTB_2_CONTROL_entropy_results$CTB_2_CONTROL_entropy_p_value[1] <-CTB_2_CONTROL_entropy_p_value
# Overwrite results file:
write.csv(CTB_2_CONTROL_entropy_results, paste0(getwd(), "/Yang_Hypernetwork_Analysis/CTB_2_CONTROL_entropy_results.csv"), row.names=FALSE)

### Plot the distribution:
hist(CTB_2_CONTROL_entropy_results$z_scores, breaks = 20, freq = FALSE, cex.lab=2, cex.axis=1.4, col = "skyblue", main = "Distribution of Z-Scores - CTB_2 CONTROL Third Trimester", xlab = "Z-Scores")
lines(density(CTB_2_CONTROL_entropy_results$z_scores), col = "red", lwd = 2)
#Add a line showing Z-score of melatonin:
abline(v =CTB_2_CONTROL_MELATONIN_entropy_Z_score, col = "green", lwd = 2)









###CTB_3
CTB_3_CONTROL_result_df <- data.frame() #Create empty data frame

for (i in 1:1000) {  # Outer loop (runs 1000 times)
  print(i)
  CTB_3_DEGs<- colnames(CTB_3_X_sd)[sample(1:4042,size=104)] # Randomly select genes from rows and label them DEGs
  cor_data_CTB_3_2 <- cor_data_CTB_3[, -na.omit(match(CTB_3_DEGs,colnames(cor_data_CTB_3)))] # This removes the random genes (DEGs) chosen in line above from the colnames
  cor_data_CTB_3_2 <- cor_data_CTB_3_2[na.omit(match(CTB_3_DEGs, rownames(cor_data_CTB_3))),] # This only keeps the random DEGs in the rownames (overall n obs/entries by 3577 columns)
  bin_CTB_3<-abs(cor_data_CTB_3_2) # Create correlation matrix w/out signs
  bin_CTB_3[which(bin_CTB_3>sd(cor_data_CTB_3_2))]<-1 #Binarise the correlation matrix
  bin_CTB_3[which(bin_CTB_3!=1)]<-0 # this creates the hypergraph incidence matrix
  hyp_CTB_3<-bin_CTB_3 %*% t(bin_CTB_3) #adjacency matrix
  CTB_3_entropy_result<- entropy(hyp_CTB_3)
  CTB_3_CONTROL_result_df <- rbind(CTB_3_CONTROL_result_df,CTB_3_entropy_result)
}
write.csv(CTB_3_CONTROL_result_df, paste0(getwd(), "/Yang_Hypernetwork_Analysis/CTB_3_CONTROL_entropy_results.csv"), row.names=FALSE)

### MELATONIN HYPERGRAPH ###

## Generates a matrix of Melatonin genes as rownames and allCTB_3 genes (- the melatonin genes of interest) as column names
CONTROL_CTB_3_DEGs<- colnames(CONTROL_X_CTB_3_sd)
CTB_3_CONTROL_melatonin_result_df <- data.frame() #Create empty data frame
cor_data_CTB_3_CONTROL <- cor_data_CTB_3[, -na.omit(match(CONTROL_CTB_3_DEGs,colnames(cor_data_CTB_3)))] #This removes the melatonin genes (DEGs) from the colnames 
cor_data_CTB_3_CONTROL <- cor_data_CTB_3_CONTROL[na.omit(match(CONTROL_CTB_3_DEGs, rownames(cor_data_CTB_3))),] # This ensures only the melatonin genes of interest are in the rows
bin_CTB_3_CONTROL<-abs(cor_data_CTB_3_CONTROL) # Create correlation matrix w/out signs
bin_CTB_3_CONTROL[which(bin_CTB_3_CONTROL>sd(cor_data_CTB_3_CONTROL))]<-1 #Binarise the correlation matrix
bin_CTB_3_CONTROL[which(bin_CTB_3_CONTROL!=1)]<-0 # this creates the hypergraph incidence matrix
hyp_CTB_3_CONTROL <-bin_CTB_3_CONTROL %*% t(bin_CTB_3_CONTROL) #adjacency matrix 
CTB_3_CONTROL_melatonin_entropy_result<- entropy(hyp_CTB_3_CONTROL)
CTB_3_CONTROL_melatonin_result_df <- rbind(CTB_3_CONTROL_melatonin_result_df,CTB_3_CONTROL_melatonin_entropy_result)

##PLOT MELATONIN HYPERGRAPH
#library(gplots)
#hm_CTB_3_CONTROL<-heatmap.2(hyp_CTB_3_CONTROL,trace="none",labRow = F,labCol=F)
#hm_CTB_3_CONTROL <- heatmap.2(hyp_CTB_3_CONTROL, trace = "none", labRow = FALSE, labCol = FALSE)
#dend_CTB_3_CONTROL<-as.hclust(hm_CTB_3_CONTROL$rowDendrogram)

### Z-scores ###
#Calculate z-scores for first hypergraphs
colnames(CTB_3_CONTROL_result_df) <- "CTB_3_CONTROL_Raw_Entropy_Values"
CTB_3_CONTROL_entropy_mean <- mean(CTB_3_CONTROL_result_df$CTB_3_CONTROL_Raw_Entropy_Values)
CTB_3_CONTROL_entropy_sd <- sd(CTB_3_CONTROL_result_df$CTB_3_CONTROL_Raw_Entropy_Values)
CTB_3_CONTROL_entropy_Z_scores <- (CTB_3_CONTROL_result_df$CTB_3_CONTROL_Raw_Entropy_Values -CTB_3_CONTROL_entropy_mean)/CTB_3_CONTROL_entropy_sd
# Add to column in entropy data frame:
z_scores_df <- data.frame(z_scores =CTB_3_CONTROL_entropy_Z_scores)
CTB_3_CONTROL_entropy_results <- cbind(CTB_3_CONTROL_result_df, z_scores_df)

### Calculate Z-score from the melatonin gene hypergraph:
CTB_3_CONTROL_MELATONIN_entropy_Z_score <- (CTB_3_CONTROL_melatonin_entropy_result-CTB_3_CONTROL_entropy_mean)/CTB_3_CONTROL_entropy_sd
### Calculate p-value and add to results file:
CTB_3_CONTROL_entropy_p_value <- 2*pnorm(-abs(CTB_3_CONTROL_MELATONIN_entropy_Z_score)) # 2* as it is two sided
#Add these to the results data frame:
CTB_3_CONTROL_entropy_results$CTB_3_CONTROL_melatonin_entropy_result <- NA
CTB_3_CONTROL_entropy_results$CTB_3_CONTROL_melatonin_entropy_result[1] <-CTB_3_CONTROL_melatonin_entropy_result
CTB_3_CONTROL_entropy_results$CONTROL_MELATONIN_entropy_Z_score <- NA
CTB_3_CONTROL_entropy_results$CONTROL_MELATONIN_entropy_Z_score[1] <-CTB_3_CONTROL_MELATONIN_entropy_Z_score
CTB_3_CONTROL_entropy_results$CTB_3_CONTROL_entropy_p_value <- NA
CTB_3_CONTROL_entropy_results$CTB_3_CONTROL_entropy_p_value[1] <-CTB_3_CONTROL_entropy_p_value
# Overwrite results file:
write.csv(CTB_3_CONTROL_entropy_results, paste0(getwd(), "/Yang_Hypernetwork_Analysis/CTB_3_CONTROL_entropy_results.csv"), row.names=FALSE)

### Plot the distribution:
hist(CTB_3_CONTROL_entropy_results$z_scores, breaks = 20, freq = FALSE, cex.lab=2, cex.axis=1.4, col = "skyblue", main = "Distribution of Z-Scores - CTB_3 CONTROL Third Trimester", xlab = "Z-Scores")
lines(density(CTB_3_CONTROL_entropy_results$z_scores), col = "red", lwd = 2)
#Add a line showing Z-score of melatonin:
abline(v =CTB_3_CONTROL_MELATONIN_entropy_Z_score, col = "green", lwd = 2)










###CTB_4
CTB_4_CONTROL_result_df <- data.frame() #Create empty data frame

for (i in 1:1000) {  # Outer loop (runs 1000 times)
  print(i)
  CTB_4_DEGs<- colnames(CTB_4_X_sd)[sample(1:4379,size=112)] # Randomly select genes from rows and label them DEGs
  cor_data_CTB_4_2 <- cor_data_CTB_4[, -na.omit(match(CTB_4_DEGs,colnames(cor_data_CTB_4)))] # This removes the random genes (DEGs) chosen in line above from the colnames
  cor_data_CTB_4_2 <- cor_data_CTB_4_2[na.omit(match(CTB_4_DEGs, rownames(cor_data_CTB_4))),] # This only keeps the random DEGs in the rownames (overall n obs/entries by 3577 columns)
  bin_CTB_4<-abs(cor_data_CTB_4_2) # Create correlation matrix w/out signs
  bin_CTB_4[which(bin_CTB_4>sd(cor_data_CTB_4_2))]<-1 #Binarise the correlation matrix
  bin_CTB_4[which(bin_CTB_4!=1)]<-0 # this creates the hypergraph incidence matrix
  hyp_CTB_4<-bin_CTB_4 %*% t(bin_CTB_4) #adjacency matrix
  CTB_4_entropy_result<- entropy(hyp_CTB_4)
  CTB_4_CONTROL_result_df <- rbind(CTB_4_CONTROL_result_df,CTB_4_entropy_result)
}
write.csv(CTB_4_CONTROL_result_df, paste0(getwd(), "/Yang_Hypernetwork_Analysis/CTB_4_CONTROL_entropy_results.csv"), row.names=FALSE)

### MELATONIN HYPERGRAPH ###

## Generates a matrix of Melatonin genes as rownames and allCTB_4 genes (- the melatonin genes of interest) as column names
CONTROL_CTB_4_DEGs<- colnames(CONTROL_X_CTB_4_sd)
CTB_4_CONTROL_melatonin_result_df <- data.frame() #Create empty data frame
cor_data_CTB_4_CONTROL <- cor_data_CTB_4[, -na.omit(match(CONTROL_CTB_4_DEGs,colnames(cor_data_CTB_4)))] #This removes the melatonin genes (DEGs) from the colnames 
cor_data_CTB_4_CONTROL <- cor_data_CTB_4_CONTROL[na.omit(match(CONTROL_CTB_4_DEGs, rownames(cor_data_CTB_4))),] # This ensures only the melatonin genes of interest are in the rows
bin_CTB_4_CONTROL<-abs(cor_data_CTB_4_CONTROL) # Create correlation matrix w/out signs
bin_CTB_4_CONTROL[which(bin_CTB_4_CONTROL>sd(cor_data_CTB_4_CONTROL))]<-1 #Binarise the correlation matrix
bin_CTB_4_CONTROL[which(bin_CTB_4_CONTROL!=1)]<-0 # this creates the hypergraph incidence matrix
hyp_CTB_4_CONTROL <-bin_CTB_4_CONTROL %*% t(bin_CTB_4_CONTROL) #adjacency matrix 
CTB_4_CONTROL_melatonin_entropy_result<- entropy(hyp_CTB_4_CONTROL)
CTB_4_CONTROL_melatonin_result_df <- rbind(CTB_4_CONTROL_melatonin_result_df,CTB_4_CONTROL_melatonin_entropy_result)

##PLOT MELATONIN HYPERGRAPH
#library(gplots)
#hm_CTB_4_CONTROL<-heatmap.2(hyp_CTB_4_CONTROL,trace="none",labRow = F,labCol=F)
#hm_CTB_4_CONTROL <- heatmap.2(hyp_CTB_4_CONTROL, trace = "none", labRow = FALSE, labCol = FALSE)
#dend_CTB_4_CONTROL<-as.hclust(hm_CTB_4_CONTROL$rowDendrogram)

### Z-scores ###
#Calculate z-scores for first hypergraphs
colnames(CTB_4_CONTROL_result_df) <- "CTB_4_CONTROL_Raw_Entropy_Values"
CTB_4_CONTROL_entropy_mean <- mean(CTB_4_CONTROL_result_df$CTB_4_CONTROL_Raw_Entropy_Values)
CTB_4_CONTROL_entropy_sd <- sd(CTB_4_CONTROL_result_df$CTB_4_CONTROL_Raw_Entropy_Values)
CTB_4_CONTROL_entropy_Z_scores <- (CTB_4_CONTROL_result_df$CTB_4_CONTROL_Raw_Entropy_Values -CTB_4_CONTROL_entropy_mean)/CTB_4_CONTROL_entropy_sd
# Add to column in entropy data frame:
z_scores_df <- data.frame(z_scores =CTB_4_CONTROL_entropy_Z_scores)
CTB_4_CONTROL_entropy_results <- cbind(CTB_4_CONTROL_result_df, z_scores_df)

### Calculate Z-score from the melatonin gene hypergraph:
CTB_4_CONTROL_MELATONIN_entropy_Z_score <- (CTB_4_CONTROL_melatonin_entropy_result-CTB_4_CONTROL_entropy_mean)/CTB_4_CONTROL_entropy_sd
### Calculate p-value and add to results file:
CTB_4_CONTROL_entropy_p_value <- 2*pnorm(-abs(CTB_4_CONTROL_MELATONIN_entropy_Z_score)) # 2* as it is two sided
#Add these to the results data frame:
CTB_4_CONTROL_entropy_results$CTB_4_CONTROL_melatonin_entropy_result <- NA
CTB_4_CONTROL_entropy_results$CTB_4_CONTROL_melatonin_entropy_result[1] <-CTB_4_CONTROL_melatonin_entropy_result
CTB_4_CONTROL_entropy_results$CONTROL_MELATONIN_entropy_Z_score <- NA
CTB_4_CONTROL_entropy_results$CONTROL_MELATONIN_entropy_Z_score[1] <-CTB_4_CONTROL_MELATONIN_entropy_Z_score
CTB_4_CONTROL_entropy_results$CTB_4_CONTROL_entropy_p_value <- NA
CTB_4_CONTROL_entropy_results$CTB_4_CONTROL_entropy_p_value[1] <-CTB_4_CONTROL_entropy_p_value
# Overwrite results file:
write.csv(CTB_4_CONTROL_entropy_results, paste0(getwd(), "/Yang_Hypernetwork_Analysis/CTB_4_CONTROL_entropy_results.csv"), row.names=FALSE)

### Plot the distribution:
hist(CTB_4_CONTROL_entropy_results$z_scores, breaks = 20, freq = FALSE, cex.lab=2, cex.axis=1.4, col = "skyblue", main = "Distribution of Z-Scores - CTB_4 CONTROL Third Trimester", xlab = "Z-Scores")
lines(density(CTB_4_CONTROL_entropy_results$z_scores), col = "red", lwd = 2)
#Add a line showing Z-score of melatonin:
abline(v =CTB_4_CONTROL_MELATONIN_entropy_Z_score, col = "green", lwd = 2)







###STB_1
STB_1_CONTROL_result_df <- data.frame() #Create empty data frame

for (i in 1:1000) {  # Outer loop (runs 1000 times)
  print(i)
  STB_1_DEGs<- colnames(STB_1_X_sd)[sample(1:4026,size=106)] # Randomly select genes from rows and label them DEGs
  cor_data_STB_1_2 <- cor_data_STB_1[, -na.omit(match(STB_1_DEGs,colnames(cor_data_STB_1)))] # This removes the random genes (DEGs) chosen in line above from the colnames
  cor_data_STB_1_2 <- cor_data_STB_1_2[na.omit(match(STB_1_DEGs, rownames(cor_data_STB_1))),] # This only keeps the random DEGs in the rownames (overall n obs/entries by 3577 columns)
  bin_STB_1<-abs(cor_data_STB_1_2) # Create correlation matrix w/out signs
  bin_STB_1[which(bin_STB_1>sd(cor_data_STB_1_2))]<-1 #Binarise the correlation matrix
  bin_STB_1[which(bin_STB_1!=1)]<-0 # this creates the hypergraph incidence matrix
  hyp_STB_1<-bin_STB_1 %*% t(bin_STB_1) #adjacency matrix
  STB_1_entropy_result<- entropy(hyp_STB_1)
  STB_1_CONTROL_result_df <- rbind(STB_1_CONTROL_result_df,STB_1_entropy_result)
}
write.csv(STB_1_CONTROL_result_df, paste0(getwd(), "/Yang_Hypernetwork_Analysis/STB_1_CONTROL_entropy_results.csv"), row.names=FALSE)

### MELATONIN HYPERGRAPH ###

## Generates a matrix of Melatonin genes as rownames and allSTB_1 genes (- the melatonin genes of interest) as column names
CONTROL_STB_1_DEGs<- colnames(CONTROL_X_STB_1_sd)
STB_1_CONTROL_melatonin_result_df <- data.frame() #Create empty data frame
cor_data_STB_1_CONTROL <- cor_data_STB_1[, -na.omit(match(CONTROL_STB_1_DEGs,colnames(cor_data_STB_1)))] #This removes the melatonin genes (DEGs) from the colnames 
cor_data_STB_1_CONTROL <- cor_data_STB_1_CONTROL[na.omit(match(CONTROL_STB_1_DEGs, rownames(cor_data_STB_1))),] # This ensures only the melatonin genes of interest are in the rows
bin_STB_1_CONTROL<-abs(cor_data_STB_1_CONTROL) # Create correlation matrix w/out signs
bin_STB_1_CONTROL[which(bin_STB_1_CONTROL>sd(cor_data_STB_1_CONTROL))]<-1 #Binarise the correlation matrix
bin_STB_1_CONTROL[which(bin_STB_1_CONTROL!=1)]<-0 # this creates the hypergraph incidence matrix
hyp_STB_1_CONTROL <-bin_STB_1_CONTROL %*% t(bin_STB_1_CONTROL) #adjacency matrix 
STB_1_CONTROL_melatonin_entropy_result<- entropy(hyp_STB_1_CONTROL)
STB_1_CONTROL_melatonin_result_df <- rbind(STB_1_CONTROL_melatonin_result_df,STB_1_CONTROL_melatonin_entropy_result)

##PLOT MELATONIN HYPERGRAPH
#library(gplots)
#hm_STB_1_CONTROL<-heatmap.2(hyp_STB_1_CONTROL,trace="none",labRow = F,labCol=F)
#hm_STB_1_CONTROL <- heatmap.2(hyp_STB_1_CONTROL, trace = "none", labRow = FALSE, labCol = FALSE)
#dend_STB_1_CONTROL<-as.hclust(hm_STB_1_CONTROL$rowDendrogram)

### Z-scores ###
#Calculate z-scores for first hypergraphs
colnames(STB_1_CONTROL_result_df) <- "STB_1_CONTROL_Raw_Entropy_Values"
STB_1_CONTROL_entropy_mean <- mean(STB_1_CONTROL_result_df$STB_1_CONTROL_Raw_Entropy_Values)
STB_1_CONTROL_entropy_sd <- sd(STB_1_CONTROL_result_df$STB_1_CONTROL_Raw_Entropy_Values)
STB_1_CONTROL_entropy_Z_scores <- (STB_1_CONTROL_result_df$STB_1_CONTROL_Raw_Entropy_Values -STB_1_CONTROL_entropy_mean)/STB_1_CONTROL_entropy_sd
# Add to column in entropy data frame:
z_scores_df <- data.frame(z_scores =STB_1_CONTROL_entropy_Z_scores)
STB_1_CONTROL_entropy_results <- cbind(STB_1_CONTROL_result_df, z_scores_df)

### Calculate Z-score from the melatonin gene hypergraph:
STB_1_CONTROL_MELATONIN_entropy_Z_score <- (STB_1_CONTROL_melatonin_entropy_result-STB_1_CONTROL_entropy_mean)/STB_1_CONTROL_entropy_sd
### Calculate p-value and add to results file:
STB_1_CONTROL_entropy_p_value <- 2*pnorm(-abs(STB_1_CONTROL_MELATONIN_entropy_Z_score)) # 2* as it is two sided
#Add these to the results data frame:
STB_1_CONTROL_entropy_results$STB_1_CONTROL_melatonin_entropy_result <- NA
STB_1_CONTROL_entropy_results$STB_1_CONTROL_melatonin_entropy_result[1] <-STB_1_CONTROL_melatonin_entropy_result
STB_1_CONTROL_entropy_results$CONTROL_MELATONIN_entropy_Z_score <- NA
STB_1_CONTROL_entropy_results$CONTROL_MELATONIN_entropy_Z_score[1] <-STB_1_CONTROL_MELATONIN_entropy_Z_score
STB_1_CONTROL_entropy_results$STB_1_CONTROL_entropy_p_value <- NA
STB_1_CONTROL_entropy_results$STB_1_CONTROL_entropy_p_value[1] <-STB_1_CONTROL_entropy_p_value
# Overwrite results file:
write.csv(STB_1_CONTROL_entropy_results, paste0(getwd(), "/Yang_Hypernetwork_Analysis/STB_1_CONTROL_entropy_results.csv"), row.names=FALSE)

### Plot the distribution:
hist(STB_1_CONTROL_entropy_results$z_scores, breaks = 20, freq = FALSE, cex.lab=2, cex.axis=1.4, col = "skyblue", main = "Distribution of Z-Scores - STB_1 CONTROL Third Trimester", xlab = "Z-Scores")
lines(density(STB_1_CONTROL_entropy_results$z_scores), col = "red", lwd = 2)
#Add a line showing Z-score of melatonin:
abline(v =STB_1_CONTROL_MELATONIN_entropy_Z_score, col = "green", lwd = 2)




