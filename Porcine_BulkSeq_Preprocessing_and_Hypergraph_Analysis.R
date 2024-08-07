
#### Analysis of Porcine Data ####

library(readr)
library(limma)
library(edgeR)
library(stringr)
library(pheatmap)
library(gplots)
library(tidyverse)
library(tidyr)
library(biomaRt)
library(dplyr)
library(umap)
library(Glimma)
library(RColorBrewer)
library(mclust)
library(rgl)
library(entropy)
library(pheatmap)

#### STEP ONE -- LOAD IN DATA AND ORGANISE ####
#Read data
setwd("C:/Users/anna-/OneDrive/Documents/Anna/MSc. Reproduction and Pregnancy/Research Project 2/pig")
data <- read.table("porcine_counts_new.txt")
summarydata <- read.table("porcine_counts_new.txt.summary")

colnames(data) <- as.character(unlist(data[1,])) #make the first row the column names
data <- data[-1,] #remove the now duplicated first row
rownames(data) <- data$Geneid # sets row names to GeneID

### Convert Ensembl IDs to gene names using BioMaRt ###
Geneid <- data$Geneid #Extract all gene names

### Using BioMaRt:
# library("biomaRt")
# listMarts() #displays all available BioMart web services
# listDatasets(ensembl) #lists all ensembl datasets
# ensembldatabase = useMart("ensembl", dataset = "hsapiens_gene_ensembl") # select hsapiens dataset
# getBM function has three arguments: filters, attributes and values.
# Filters restricts the query (use listFilters(ensembldatabase))
# Attributes define the values you want to retrieve/the output (use listAttributes(ensembldatabase))
# Values are the list of genes it is searching with
# getLDS (Get linked dataset) is used to link two biomart datasets
# AttributesL and filtersL are used to specify features of the linked dataset.

###Connect to porcine Ensembl dataset (use archived host to avoid server error):
ensembldatabase <- useMart(biomart="ensembl", dataset = "sscrofa_gene_ensembl", verbose=TRUE, host="https://feb2021.archive.ensembl.org")

#Set up a query to convert the gene names using the get BM function (use getBM when only using one database, get LDS for two):
gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters ="ensembl_gene_id",
                    values = Geneid,
                    mart = ensembldatabase)

#GENENAME = external name
#GENEID = ensembl ID

colnames(gene_names)[1] <- 'GENEID'
colnames(gene_names)[2] <- 'GENENAME'

gene_names$GENENAME <-as.character(gene_names$GENENAME)
gene_names$GENENAME <- make.unique(gene_names$GENENAME,sep = "~")

data <- data[row.names(data) %in% gene_names$GENEID,]#remove any genes that don't match with identifier 
genes <- gene_names$GENENAME# takes first column as gene name

# Check if data$Geneid is in the same order as gene_names$GENEID
is_same_order <- all(identical(data$Geneid, gene_names$GENEID))
# Print the result
print(is_same_order)

#This was false, so need to use match function:
gene_names <- gene_names[match(data$Geneid, gene_names$GENEID), ]
#check again:
is_same_order <- all(identical(data$Geneid, gene_names$GENEID))
print(is_same_order)

#merge the data frames based on Geneid column
merged_data <- merge(data, gene_names, by.x = "Geneid", by.y = "GENEID", all.x = TRUE)
#Move GENENAME column to the beginning
merged_data <- cbind(merged_data[ncol(merged_data)], merged_data[-ncol(merged_data)])

data1 <- merged_data

genes <- data1$GENENAME# takes first column as gene name
rownames(data1) <- make.unique(genes,sep = "~") #reimput geneID as row names. This will add .1 to the end of gene names

drop <- c("Chr", "Start", "End", "Strand", "Length", "Geneid") # removing stuff you don't need
data1 <- data1[, !(names(data1) %in% drop)]

#rename columns to just sample names (keep .sorted.bam columns only):
sample_names <- strsplit(colnames(data1), "/")
sample_names <- unlist(lapply(sample_names[-1], "[[", 9))# select the element you want ignoring the first columns
#keep only first bit before full stop:
sample_names <- sub("\\..*", "", sample_names)

#set col names of data
colnames(data1) <- c("GENENAME", sample_names)

remove_rows_with_tilde <- function(df) {
  df[!grepl("~", rownames(df)), ]}

data1 <- remove_rows_with_tilde(data1)
data1 <- data.frame(data1, row.names = 1)

#Check if the colname are there and assigns them to melatonin, if not then they will be assigned to control
sample_info <- read.csv("porcine_sample_info.csv")

#need to reorder count data cols to match the sample_info data or other way around
data1 <- data1[,c(1, na.omit(match(sample_info$run,colnames(data1))))]

#delete second row
data1 <- data1[,-2]

## Overall data should have sample number as colnames and genes as rownames ###



##### STEP TWO -- CREATE YOUR COMPARISON GROUPS #####
# Using edgeR:
# The code below uses the glm approach
# There are two testing methods: Liklihood ratio and quasi-likelihood (QL)
# The QL method is recommended for DE analysis on bulk seq
# The likelihood ratio test can be used in scRNA seq (and data with no replicates)

# The group needs to be a list of values e.g. group<-factor(c(1,1,2,2))

# specify which column describes your groups #
group <- sample_info$group
design <- model.matrix(~0+group) #This is a design matrix without an intercept (~0). 

### CHECK THAT THESE ARE THE CORRECT WAY ROUND --- = case - control
# specify name of group comparisons and the column numbers to compare #all standard EdgeR stuff
#make contrasts generates a contrast matrix for the coefficients in a linear model
#This one is comparing between the melatonin and control groups
my_contrast <- makeContrasts(melatoninvscontrol = groupmelatonin - groupcontrol,levels= design) 

######CHECK THIS IS THE RIGHT WAY AROUND!!!
groups <- group %in% c("control","melatonin") #control goes first in list, change according to your data and change for each analysis


##### STEP THREE -- FILTER AND NORMALISE DATA #####
data1[, 1] <- as.numeric(data1[, 1]) #converts data to numeric form
data1[, 2] <- as.numeric(data1[, 2])
data1[, 3] <- as.numeric(data1[, 3])
data1[, 4] <- as.numeric(data1[, 4])

edgeRformat_data <- DGEList(counts=data1, group=group) # DGElist is the data object that is being analysed with a grouping factor
keep <- filterByExpr(edgeRformat_data) #filters out lowly expressed
edgeRformat_data <- edgeRformat_data[keep,,keep.lib.sizes=F] #more filtering
edgeRformat_data <- calcNormFactors(edgeRformat_data) #normalises for library sizes
edgeRformat_data <- estimateDisp(edgeRformat_data, design) #recommended to estimate common dispersion and tagwise in one go
fit <- glmQLFit(edgeRformat_data, design) # this fits the negative binomial GLM for each tag(aka DEG) and produces a new object DGEGLM


##### STEP FOUR -- RUN THE DIFFERENTIAL GENE EXPRESSION FOR YOUR COMPARISON GROUP #####
# run this for different groups #

qlf.2v1 <- glmQLFTest(fit, contrast = my_contrast[,"melatoninvscontrol"]) # dispersion estimation and hypothesis testing (the QL Test)
dge_output <- topTags(qlf.2v1, n=Inf, p.value = 0.05)$table

#### CONVERT PORCINE DEGs TO HUMAN ####
#Extract DEG names from the rownames
porcine_gene_names <- rownames(dge_output)

#Connect to Ensembl datasets (use old archive to avoid server error)
human <- useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose=TRUE, host="https://feb2021.archive.ensembl.org")
porcine <- useMart(biomart="ensembl", dataset = "sscrofa_gene_ensembl", verbose=TRUE, host="https://feb2021.archive.ensembl.org")

#Set up a query to convert the gene names using the get LDS function:
genesV2 <- getLDS(attributes = c("external_gene_name"),
                  filters ="external_gene_name",
                  values = porcine_gene_names,
                  mart = porcine,
                  attributesL = c("hgnc_symbol"),
                  martL = human,
                  uniqueRows = TRUE)

#Check for duplicates:
duplicates <- genesV2[duplicated(genesV2$HGNC.symbol), ]
duplicates

# Match row names (GeneID) in dge_output with GeneID in genesV2
matched_genes <- intersect(rownames(dge_output), genesV2$Gene.name) #These are the pig genes that have a human ortholog

# Add HGNC.symbol column to dge_output based on matching row names
dge_output$HGNC.symbol <- NA  # make a new column
match_indices <- match(rownames(dge_output), genesV2$Gene.name)
dge_output$HGNC.symbol <- genesV2$HGNC.symbol[match_indices]
dge_output<- dge_output[!is.na(dge_output$HGNC.symbol), ]#omit rows without human ortholog

#write.csv(dge_output, row.names=TRUE, "C:/Users/anna-/OneDrive/Documents/Anna/MSc. Reproduction and Pregnancy/Research Project 2/pig/DEGs_MELvsNormal.csv")
#write.csv(data1, row.names=TRUE, "~/RP2/pig/new/processedporcinedata.csv")




#### PCA ####

# Get log2 counts per million
logcounts <- cpm(edgeRformat_data,log = T) #cpm is counts per million. this is edgeR function

#write.csv(logcounts, "porcine_logcounts.csv")

## ONLY NEED TO RUN ONCE --- PCA for all of the data - not for the analysis chosen ##
pcadata_all<-pca(t(logcounts),ncomp = 2) #pca analysis, ncomp is number of principal components

#calculate variance
eigs <- pcadata_all$sdev^2 #calculates the variance (eigenvalues) of the PCs
PCA1_all_var <- eigs[1] / sum(eigs) #calculates the proportion of total variance explained by first PC
PCA1_all_var<-round(PCA1_all_var,digits = 2)*100 #rounds to 2 decimal places and converts to %

eigs <- pcadata_all$sdev^2 
PCA2_all_var <- eigs[2] / sum(eigs) #calculates the proportion of total variance explained by second PC
PCA2_all_var<-round(PCA2_all_var,digits = 2)*100 #rounds and converts to %

pcadata_plot <- data.frame(pcadata_all$variates$X[,1:2], 
                           group=sample_info$group) #extracts pca scores, adds column called group

#plotting PCs and saving to device
#svg("C:/Users/anna-/OneDrive/Documents/Anna/MSc. Reproduction and Pregnancy/Research Project 2/pig/pig_analysis_PCA.svg", height = 10, width = 10)
#ggplot(data=pcadata_plot,aes(x=PC1,y=PC2,col=group,))+
# geom_point(size=5)+
# theme_bw(base_size = 22)+
# labs(title= "PCA for all samples")+
# xlab(label= paste("PCA1 (variance = ", PCA1_all_var, "%)",sep=""))+
# ylab(label= paste("PCA2 (variance = ", PCA2_all_var, "%)",sep=""))
#dev.off() #this closes the plot window


#plot MDS using Glimma - this will open a web browser page and produce an interactive graph
#MDS (multidimensional scaling) is used to map data into 2D space
#labels <- paste(sample_info$group, sample_info$sample_name)
#group <- sample_info$group
#glMDSPlot(logcounts, labels=labels, groups=group, folder="mds")




#### Hypergraph: Porcine Transcriptome vs. DEGs ####

data1<- read.csv("~/RP2/pig/new/processedporcinedata.csv")

#First ensure that data1 is in matrix form:
data1 <- as.matrix(data1)

#check genes are col headings, samples are row headings
#data1 <- t(data1)

#Remove SD's of 0 and rename df
sd.scores1<-apply(data1,2,sd)
data2<-data1[,which(sd.scores1>0)]

#Correlation
cor_data_allporcine <- cor(data2, data2)

## Construct the Hypergraph:
#As there are only 4 samples, it would be inappropriate to use +/- 1SD to binarise the matrix
#In cases of low n, use a percentile cut off. Aim to get a few hundred genes. Should be between 90 and 95%.

PORCINE_DEGs <- genesV2$Gene.name #character vector of DEGs
result_df <- data.frame() #Create empty data frame
cor_data_allvsdegs <- cor_data_allporcine[,-na.omit(match(PORCINE_DEGs, colnames(cor_data_allporcine)))] #This removes the DEGs from the colnames of correlation matrix
cor_data_allvsdegs <- cor_data_allvsdegs[na.omit(match(PORCINE_DEGs, rownames(cor_data_allporcine))),] # This ensures only the DEGs are in the rows
bin_allvsdegs<-abs(cor_data_allvsdegs) # Take absolute values

bin_allvsdegs[which(bin_allvsdegs >= quantile(cor_data_allvsdegs, 0.95))]<-1 #Binarize based on percentile of r-value
bin_allvsdegs[which(bin_allvsdegs!=1)]<-0 # all other values assume 0 - this is the incidence matrix

hyp_allvsdegs <-bin_allvsdegs %*% t(bin_allvsdegs) #multiply by transpose - this is the adjacency matrix 
allvsdegs_entropy_result<- entropy(hyp_allvsdegs)
result_df <- rbind(result_df,allvsdegs_entropy_result)

## PLOT THE HEATMAP ##
allvsdegs_hm <- pheatmap::pheatmap(hyp_allvsdegs, fontsize_row = 2,fontsize_col = 2,cutree_rows = 2,cutree_cols = 2)
allvsdegs_hm_clusters <- as.data.frame(cutree(allvsdegs_hm$tree_row,2))
#optional creating column with the names of the nodes
allvsdegs_hm_clusters$Gene <- rownames(allvsdegs_hm_clusters)



## CORRESPONDENCE ANALYSIS ##
#Separating out into 2 dataframes, each containing names of nodes in that cluster
#I just check the number of datapoints to see which cluster is which (normally pretty obvious from heatmap)
allvsdegs_hm_cluster_1 <- as.data.frame(allvsdegs_hm_clusters[which(allvsdegs_hm_clusters$`cutree(allvsdegs_hm$tree_row, 2)` == 2),])
allvsdegs_hm_cluster_2 <- as.data.frame(allvsdegs_hm_clusters[which(allvsdegs_hm_clusters$`cutree(allvsdegs_hm$tree_row, 2)` == 1),])
## Overall - this is a list of the genes in the central cluster
# Put this into a character vector
allvsdegs_central_cluster<- allvsdegs_hm_cluster_2$Gene
#Restrict the rows to only include genes in central cluster
row_restricted_allvsdegs <- bin_allvsdegs[na.omit(match(allvsdegs_central_cluster,rownames(bin_allvsdegs))),]
#Restrict columns to only include genes from STB_1 that are correlated with every central cluster gene
galois_allvsdegs <- row_restricted_allvsdegs [,which(colSums(row_restricted_allvsdegs)==length(allvsdegs_central_cluster))]
galois_allvsdegs <- as.data.frame(galois_allvsdegs)
galois_allvsdegs <- tibble::rownames_to_column(galois_allvsdegs)



#### expression of key genes (logcounts) ####
#convert logcounts from matrix to DataFrame:
logcounts <- as.data.frame(logcounts)

genes_of_interest <- c("MTNR1A", "MTNR1B", "AANAT", "ASMT")
subset_logcounts <- logcounts[rownames(logcounts) %in% genes_of_interest, ]
print(subset_logcounts)
