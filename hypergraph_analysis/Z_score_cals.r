            
          ### Calculate Z-scores ###

## Load in packages:
library(dplyr)

setwd("~/RP2/vento-tormo/VT_all_anndata")

# Load in all data

folder_path <- paste0(getwd(), "Entropy_Results/")
all_files <- list.files(folder_path, pattern="\\.csv$", full.names=TRUE)
# Loop through each file and assign to the global environment
for (file in files) {
  # Create a variable name from the file name (without extension)
  var_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the file
  data <- read.csv(file, header = TRUE)
  
  # Assign the dataframe to a variable with the file name
  assign(var_name, data, envir = .GlobalEnv)
}


# Change column names
colnames(CTB_1_MTNR1A_entropy_results) <- "CTB_1_MTNR1A_Raw_Entropy_Values"
colnames(CTB_2_MTNR1A_entropy_results) <- "CTB_2_MTNR1A_Raw_Entropy_Values"
colnames(CTB_3_MTNR1A_entropy_results) <- "CTB_3_MTNR1A_Raw_Entropy_Values"
colnames(CTB_4_MTNR1A_entropy_results) <- "CTB_4_MTNR1A_Raw_Entropy_Values"
colnames(END_1_MTNR1A_entropy_results) <- "END_1_MTNR1A_Raw_Entropy_Values"
colnames(EVT_1_MTNR1A_entropy_results) <- "EVT_1_MTNR1A_Raw_Entropy_Values"
colnames(EVT_2_MTNR1A_entropy_results) <- "EVT_2_MTNR1A_Raw_Entropy_Values"
colnames(EVT_3_MTNR1A_entropy_results) <- "EVT_3_MTNR1A_Raw_Entropy_Values"
colnames(FB_1_MTNR1A_entropy_results) <- "FB_1_MTNR1A_Raw_Entropy_Values"
colnames(FB_2_MTNR1A_entropy_results) <- "FB_2_MTNR1A_Raw_Entropy_Values"
colnames(FB_3_MTNR1A_entropy_results) <- "FB_3_MTNR1A_Raw_Entropy_Values"
colnames(HB_1_MTNR1A_entropy_results) <- "HB_1_MTNR1A_Raw_Entropy_Values"
colnames(HB_2_MTNR1A_entropy_results) <- "HB_2_MTNR1A_Raw_Entropy_Values"
colnames(HB_3_MTNR1A_entropy_results) <- "HB_3_MTNR1A_Raw_Entropy_Values"
colnames(MIC_1_MTNR1A_entropy_results) <- "MIC_1_MTNR1A_Raw_Entropy_Values"
colnames(MIC_2_MTNR1A_entropy_results) <- "MIC_2_MTNR1A_Raw_Entropy_Values"
colnames(MIC_3_MTNR1A_entropy_results) <- "MIC_3_MTNR1A_Raw_Entropy_Values"
colnames(SM_1_MTNR1A_entropy_results) <- "SM_1_MTNR1A_Raw_Entropy_Values"
colnames(SM_2_MTNR1A_entropy_results) <- "SM_2_MTNR1A_Raw_Entropy_Values"
colnames(STB_1_MTNR1A_entropy_results) <- "STB_1_MTNR1A_Raw_Entropy_Values"


# Calculate mean and SD
CTB_1_MTNR1A_entropy_mean <- mean(CTB_1_MTNR1A__entropy_results$CTB_1_MTNR1A_Raw_Entropy_Values)
CTB_1_MTNR1A_entropy_sd <- sd(CTB_1_MTNR1A_entropy_results$CTB_1_MTNR1A_Raw_Entropy_Values)

CTB_2_MTNR1A_entropy_mean <- mean(CTB_2_MTNR1A_entropy_results$CTB_2_MTNR1A_Raw_Entropy_Values)
CTB_2_MTNR1A_entropy_sd <- sd(CTB_2_MTNR1A_entropy_results$CTB_2_MTNR1A_Raw_Entropy_Values)

CTB_3_MTNR1A_entropy_mean <- mean(CTB_3_MTNR1A_entropy_results$CTB_3_MTNR1A_Raw_Entropy_Values)
CTB_3_MTNR1A_entropy_sd <- sd(CTB_3_MTNR1A_entropy_results$CTB_3_MTNR1A_Raw_Entropy_Values)

CTB_4_MTNR1A_entropy_mean <- mean(CTB_4_MTNR1A_entropy_results$CTB_4_MTNR1A_Raw_Entropy_Values)
CTB_4_MTNR1A_entropy_sd <- sd(CTB_4_MTNR1A_entropy_results$CTB_4_MTNR1A_Raw_Entropy_Values)

END_1_MTNR1A_entropy_mean <- mean(END_1_MTNR1A_entropy_results$END_1_MTNR1A_Raw_Entropy_Values)
END_1_MTNR1A_entropy_sd <- sd(END_1_MTNR1A_entropy_results$END_1_MTNR1A_Raw_Entropy_Values)

EVT_1_MTNR1A_entropy_mean <- mean(EVT_1_MTNR1A_entropy_results$EVT_1_MTNR1A_Raw_Entropy_Values)
EVT_1_MTNR1A_entropy_sd <- sd(EVT_1_MTNR1A_entropy_results$EVT_1_MTNR1A_Raw_Entropy_Values)

EVT_2_MTNR1A_entropy_mean <- mean(EVT_2_MTNR1A_entropy_results$EVT_2_MTNR1A_Raw_Entropy_Values)
EVT_2_MTNR1A_entropy_sd <- sd(EVT_2_MTNR1A_entropy_results$EVT_2_MTNR1A_Raw_Entropy_Values)

EVT_3_MTNR1A_entropy_mean <- mean(EVT_3_MTNR1A_entropy_results$EVT_3_MTNR1A_Raw_Entropy_Values)
EVT_3_MTNR1A_entropy_sd <- sd(EVT_3_MTNR1A_entropy_results$EVT_3_MTNR1A_Raw_Entropy_Values)

FB_1_MTNR1A_entropy_mean <- mean(FB_1_MTNR1A_entropy_results$FB_1_MTNR1A_Raw_Entropy_Values)
FB_1_MTNR1A_entropy_sd <- sd(FB_1_MTNR1A_entropy_results$FB_1_MTNR1A_Raw_Entropy_Values)

FB_2_MTNR1A_entropy_mean <- mean(FB_2_MTNR1A_entropy_results$FB_2_MTNR1A_Raw_Entropy_Values)
FB_2_MTNR1A_entropy_sd <- sd(FB_2_MTNR1A_entropy_results$FB_2_MTNR1A_Raw_Entropy_Values)

FB_3_MTNR1A_entropy_mean <- mean(FB_3_MTNR1A_entropy_results$FB_3_MTNR1A_Raw_Entropy_Values)
FB_3_MTNR1A_entropy_sd <- sd(FB_3_MTNR1A_entropy_results$FB_3_MTNR1A_Raw_Entropy_Values)

HB_1_MTNR1A_entropy_mean <- mean(HB_1_MTNR1A_entropy_results$HB_1_MTNR1A_Raw_Entropy_Values)
HB_1_MTNR1A_entropy_sd <- sd(HB_1_MTNR1A_entropy_results$HB_1_MTNR1A_Raw_Entropy_Values)

HB_2_MTNR1A_entropy_mean <- mean(HB_2_MTNR1A_entropy_results$HB_2_MTNR1A_Raw_Entropy_Values)
HB_2_MTNR1A_entropy_sd <- sd(HB_2_MTNR1A_entropy_results$HB_2_MTNR1A_Raw_Entropy_Values)

HB_3_MTNR1A_entropy_mean <- mean(HB_3_MTNR1A_entropy_results$HB_3_MTNR1A_Raw_Entropy_Values)
HB_3_MTNR1A_entropy_sd <- sd(HB_3_MTNR1A_entropy_results$HB_3_MTNR1A_Raw_Entropy_Values)

MIC_1_MTNR1A_entropy_mean <- mean(MIC_1_MTNR1A_entropy_results$MIC_1_MTNR1A_Raw_Entropy_Values)
MIC_1_MTNR1A_entropy_sd <- sd(MIC_1_MTNR1A_entropy_results$MIC_1_MTNR1A_Raw_Entropy_Values)

MIC_2_MTNR1A_entropy_mean <- mean(MIC_2_MTNR1A_entropy_results$MIC_2_MTNR1A_Raw_Entropy_Values)
MIC_2_MTNR1A_entropy_sd <- sd(MIC_2_MTNR1A_entropy_results$MIC_2_MTNR1A_Raw_Entropy_Values)

MIC_3_MTNR1A_entropy_mean <- mean(MIC_3_MTNR1A_entropy_results$MIC_3_MTNR1A_Raw_Entropy_Values)
MIC_3_MTNR1A_entropy_sd <- sd(MIC_3_MTNR1A_entropy_results$MIC_3_MTNR1A_Raw_Entropy_Values)

SM_1_MTNR1A_entropy_mean <- mean(SM_1_MTNR1A_entropy_results$SM_1_MTNR1A_Raw_Entropy_Values)
SM_1_MTNR1A_entropy_sd <- sd(SM_1_MTNR1A_entropy_results$SM_1_MTNR1A_Raw_Entropy_Values)

SM_2_MTNR1A_entropy_mean <- mean(SM_2_MTNR1A_entropy_results$SM_2_MTNR1A_Raw_Entropy_Values)
SM_2_MTNR1A_entropy_sd <- sd(SM_2_MTNR1A_entropy_results$SM_2_MTNR1A_Raw_Entropy_Values)

STB_1_MTNR1A_entropy_mean <- mean(STB_1_MTNR1A_entropy_results$STB_1_MTNR1A_Raw_Entropy_Values)
STB_1_MTNR1A_entropy_sd <- sd(STB_1_MTNR1A_entropy_results$STB_1_MTNR1A_Raw_Entropy_Values)

                ### CALCULATE Z-SCORES ###
### CTB_1:
# Calculate Z scores
CTB_1_MTNR1Aentropy_Z_scores <- (CTB_1_MTNR1A_entropy_results$CTB_1_MTNR1A_Raw_Entropy_Values - CTB_1_MTNR1A_entropy_mean)/CTB_1_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = CTB_1_MTNR1A_entropy_Z_scores)
CTB_1_MTNR1A_entropy_results <- cbind(CTB_1_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(CTB_1_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/CTB_1_MTNR1A_entropy_results.csv"), row.names=FALSE)

### CTB_2:
# Calculate Z scores
CTB_2_MTNR1A_entropy_Z_scores <- (CTB_2_MTNR1A_entropy_results$CTB_2_MTNR1A_Raw_Entropy_Values - CTB_2_MTNR1A_entropy_mean)/CTB_2_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = CTB_2_MTNR1A_entropy_Z_scores)
CTB_2_MTNR1A_entropy_results <- cbind(CTB_2_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(CTB_2_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/CTB_2_MTNR1A_entropy_results.csv"), row.names=FALSE)

### CTB_3:
# Calculate Z scores
CTB_3_MTNR1A_entropy_Z_scores <- (CTB_3_MTNR1A_entropy_results$CTB_3_MTNR1A_Raw_Entropy_Values - CTB_3_MTNR1A_entropy_mean)/CTB_3_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = CTB_3_MTNR1A_entropy_Z_scores)
CTB_3_MTNR1A_entropy_results <- cbind(CTB_3_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(CTB_3_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/CTB_3_MTNR1A_entropy_results.csv"), row.names=FALSE)

### CTB_4:
# Calculate Z scores
CTB_4_MTNR1A_entropy_Z_scores <- (CTB_4_MTNR1A_entropy_results$CTB_4_MTNR1A_Raw_Entropy_Values - CTB_4_MTNR1A_entropy_mean)/CTB_4_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = CTB_4_MTNR1A_entropy_Z_scores)
CTB_4_MTNR1A_entropy_results <- cbind(CTB_4_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(CTB_4_MTNR1A_entropy_results,paste0(getwd(), "/Entropy_Results/CTB_4_MTNR1A_entropy_results.csv"), row.names=FALSE)

### END_1:
# Calculate Z scores
END_1_MTNR1A_entropy_Z_scores <- (END_1_MTNR1A_entropy_results$END_1_MTNR1A_Raw_Entropy_Values - END_1_MTNR1A_entropy_mean)/END_1_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = END_1_MTNR1A_entropy_Z_scores)
END_1_MTNR1A_entropy_results <- cbind(END_1_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(END_1_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/END_1_MTNR1A_entropy_results.csv"), row.names=FALSE)

### EVT_1:
# Calculate Z scores
EVT_1_MTNR1A_entropy_Z_scores <- (EVT_1_MTNR1A_entropy_results$EVT_1_MTNR1A_Raw_Entropy_Values - EVT_1_MTNR1A_entropy_mean)/EVT_1_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = EVT_1_MTNR1A_entropy_Z_scores)
EVT_1_MTNR1A_entropy_results <- cbind(EVT_1_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(EVT_1_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/EVT_1_MTNR1A_entropy_results.csv"), row.names=FALSE)

### EVT_2:
# Calculate Z scores
EVT_2_MTNR1A_entropy_Z_scores <- (EVT_2_MTNR1A_entropy_results$EVT_2_MTNR1A_Raw_Entropy_Values - EVT_2_MTNR1A_entropy_mean)/EVT_2_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = EVT_2_MTNR1A_entropy_Z_scores)
EVT_2_MTNR1A_entropy_results <- cbind(EVT_2_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(EVT_2_MTNR1A_entropy_results,paste0(getwd(), "/Entropy_Results/EVT_2_MTNR1A_entropy_results.csv"), row.names=FALSE)

### EVT_3:
# Calculate Z scores
EVT_3_MTNR1A_entropy_Z_scores <- (EVT_3_MTNR1A_entropy_results$EVT_3_MTNR1A_Raw_Entropy_Values - EVT_3_MTNR1A_entropy_mean)/EVT_3_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = EVT_3_MTNR1A_entropy_Z_scores)
EVT_3_MTNR1A_entropy_results <- cbind(EVT_3_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(EVT_3_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/EVT_3_MTNR1A_entropy_results.csv"), row.names=FALSE)

### FB_1:
# Calculate Z scores
FB_1_MTNR1A_entropy_Z_scores <- (FB_1_MTNR1A_entropy_results$FB_1_MTNR1A_Raw_Entropy_Values - FB_1_MTNR1A_entropy_mean)/FB_1_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = FB_1_MTNR1A_entropy_Z_scores)
FB_1_MTNR1A_entropy_results <- cbind(FB_1_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(FB_1_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/FB_1_MTNR1A_entropy_results.csv"), row.names=FALSE)

### FB_2:
# Calculate Z scores
FB_2_MTNR1A_entropy_Z_scores <- (FB_2_MTNR1A_entropy_results$FB_2_MTNR1A_Raw_Entropy_Values - FB_2_MTNR1A_entropy_mean)/FB_2_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = FB_2_MTNR1A_entropy_Z_scores)
FB_2_MTNR1A_entropy_results <- cbind(FB_2_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(FB_2_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/FB_2_MTNR1A_entropy_results.csv"), row.names=FALSE)

### FB_3:
# Calculate Z scores
FB_3_MTNR1A_entropy_Z_scores <- (FB_3_MTNR1A_entropy_results$FB_3_MTNR1A_Raw_Entropy_Values - FB_3_MTNR1A_entropy_mean)/FB_3_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = FB_3_MTNR1A_entropy_Z_scores)
FB_3_MTNR1A_entropy_results <- cbind(FB_3_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(FB_3_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/FB_3_MTNR1A_entropy_results.csv"), row.names=FALSE)

### HB_1:
# Calculate Z scores
HB_1_MTNR1A_entropy_Z_scores <- (HB_1_MTNR1A_entropy_results$HB_1_MTNR1A_Raw_Entropy_Values - HB_1_MTNR1A_entropy_mean)/HB_1_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = HB_1_MTNR1A_entropy_Z_scores)
HB_1_MTNR1A_entropy_results <- cbind(HB_1_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(HB_1_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/HB_1_MTNR1A_entropy_results.csv"), row.names=FALSE)

### HB_2:
# Calculate Z scores
HB_2_MTNR1A_entropy_Z_scores <- (HB_2_MTNR1A_entropy_results$HB_2_MTNR1A_Raw_Entropy_Values - HB_2_MTNR1A_entropy_mean)/HB_2_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = HB_2_MTNR1A_entropy_Z_scores)
HB_2_MTNR1A_entropy_results <- cbind(HB_2_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(HB_2_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/HB_2_MTNR1A_entropy_results.csv"), row.names=FALSE)

### HB_3:
# Calculate Z scores
HB_3_MTNR1A_entropy_Z_scores <- (HB_3_MTNR1A_entropy_results$HB_3_MTNR1A_Raw_Entropy_Values - HB_3_MTNR1A_entropy_mean)/HB_3_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = HB_3_MTNR1A_entropy_Z_scores)
HB_3_MTNR1A_entropy_results <- cbind(HB_3_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(HB_3_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/HB_3_MTNR1A_entropy_results.csv"), row.names=FALSE)

### MIC_1:
# Calculate Z scores
MIC_1_MTNR1A_entropy_Z_scores <- (MIC_1_MTNR1A_entropy_results$MIC_1_MTNR1A_Raw_Entropy_Values - MIC_1_MTNR1A_entropy_mean)/MIC_1_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = MIC_1_MTNR1A_entropy_Z_scores)
MIC_1_MTNR1A_entropy_results <- cbind(MIC_1_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(MIC_1_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/MIC_1_MTNR1A_entropy_results.csv"), row.names=FALSE)

### MIC_2:
# Calculate Z scores
MIC_2_MTNR1A_entropy_Z_scores <- (MIC_2_MTNR1A_entropy_results$MIC_2_MTNR1A_Raw_Entropy_Values - MIC_2_MTNR1A_entropy_mean)/MIC_2_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = MIC_2_MTNR1A_entropy_Z_scores)
MIC_2_MTNR1A_entropy_results <- cbind(MIC_2_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(MIC_2_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/MIC_2_MTNR1A_entropy_results.csv"), row.names=FALSE)

### MIC_3:
# Calculate Z scores
MIC_3_MTNR1A_entropy_Z_scores <- (MIC_3_MTNR1A_entropy_results$MIC_3_MTNR1A_Raw_Entropy_Values - MIC_3_MTNR1A_entropy_mean)/MIC_3_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = MIC_3_MTNR1A_entropy_Z_scores)
MIC_3_MTNR1A_entropy_results <- cbind(MIC_3_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(MIC_3_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/MIC_3_MTNR1A_entropy_results.csv"), row.names=FALSE)

### SM_1:
# Calculate Z scores
SM_1_MTNR1A_entropy_Z_scores <- (SM_1_MTNR1A_entropy_results$SM_1_MTNR1A_Raw_Entropy_Values - SM_1_MTNR1A_entropy_mean)/SM_1_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = SM_1_MTNR1A_entropy_Z_scores)
SM_1_MTNR1A_entropy_results <- cbind(SM_1_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(SM_1_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/SM_1_MTNR1A_entropy_results.csv"), row.names=FALSE)

### SM_2:
# Calculate Z scores
SM_2_MTNR1A_entropy_Z_scores <- (SM_2_MTNR1A_entropy_results$SM_2_MTNR1A_Raw_Entropy_Values - SM_2_MTNR1A_entropy_mean)/SM_2_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = SM_2_MTNR1A_entropy_Z_scores)
SM_2_MTNR1A_entropy_results <- cbind(SM_2_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(SM_2_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/SM_2_MTNR1A_entropy_results.csv"), row.names=FALSE)

### STB_1:
# Calculate Z scores
STB_1_MTNR1A_entropy_Z_scores <- (STB_1_MTNR1A_entropy_results$STB_1_MTNR1A_Raw_Entropy_Values - STB_1_MTNR1A_entropy_mean)/STB_1_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = STB_1_MTNR1A_entropy_Z_scores)
STB_1_MTNR1A_entropy_results <- cbind(STB_1_MTNR1A_entropy_results, z_scores_df)
# Overwrite:
write.csv(STB_1_MTNR1A_entropy_results, paste0(getwd(), "/Entropy_Results/STB_1_MTNR1A_entropy_results.csv"), row.names=FALSE)


# Merge all z-scores into one CSV:
Combined_MTNR1A_z_scores <-  cbind(CTB_1_MTNR1A = CTB_1_MTNR1A_entropy_Z_scores, CTB_2_MTNR1A = CTB_2_MTNR1A_entropy_Z_scores, CTB_3_MTNR1A = CTB_3_MTNR1A_entropy_Z_scores,
                                 CTB_4_MTNR1A = CTB_4_MTNR1A_entropy_Z_scores, END_1_MTNR1A = END_1_MTNR1A_entropy_Z_scores, EVT_1_MTNR1A = EVT_1_MTNR1A_entropy_Z_scores,
                                 EVT_2_MTNR1A = EVT_2_MTNR1A_entropy_Z_scores, EVT_3_MTNR1A = EVT_3_MTNR1A_entropy_Z_scores, FB_1_MTNR1A = FB_1_MTNR1A_entropy_Z_scores,
                                 FB_2_MTNR1A = FB_2_MTNR1A_entropy_Z_scores, FB_3_MTNR1A = FB_3_MTNR1A_entropy_Z_scores, HB_1_MTNR1A = HB_1_MTNR1A_entropy_Z_scores,
                                 HB_2_MTNR1A = HB_2_MTNR1A_entropy_Z_scores, HB_3_MTNR1A = HB_3_MTNR1A_entropy_Z_scores, MIC_1_MTNR1A = MIC_1_MTNR1A_entropy_Z_scores,
                                 MIC_2_MTNR1A = MIC_2_MTNR1A_entropy_Z_scores, MIC_3_MTNR1A = MIC_3_MTNR1A_entropy_Z_scores, SM_1_MTNR1A = SM_1_MTNR1A_entropy_Z_scores,SM_2_MTNR1A = SM_2_MTNR1A_entropy_Z_scores,
                                 STB_1_MTNR1A = STB_1_MTNR1A_entropy_Z_scores)

write.csv(Combined_MTNR1A_z_scores, paste0(getwd(), "/Entropy_Results/Combined_MTNR1A_z_scores.csv"), row.names=FALSE)



            ### GROUPING ALL DATA BY CELL TYPE ###

library(dplyr)
library(tidyr)

ALL_CTB_MTNR1A_Raw_entropy <- bind_cols(CTB_1_MTNR1A_entropy_results$CTB_1_Raw_Entropy_Values, CTB_2_MTNR1A_entropy_results$CTB_2_Raw_Entropy_Values, CTB_3_MTNR1A_entropy_results$CTB_3_Raw_Entropy_Values, CTB_4_MTNR1A_entropy_results$CTB_4_Raw_Entropy_Values)
# Reshape the combined data into a single column
ALL_CTB_MTNR1A_entropy <- ALL_CTB_MTNR1A_entropy %>%
  pivot_longer(everything(), names_to = "CTB_Group", values_to = "CTB_Entropy")

ALL_END_MTNR1A_Raw_entropy <- bind_cols(END_1_MTNR1A_entropy_results$END_1_Raw_Entropy_Values)
# Reshape the combined data into a single column
ALL_END_MTNR1A_Raw_entropy <- ALL_END_Raw_MTNR1A_entropy %>%
  pivot_longer(everything(), names_to = "END_Group", values_to = "END_Entropy")

ALL_EVT_MTNR1A_Raw_entropy <- bind_cols(EVT_1_MTNR1A_entropy_results$EVT_1_Raw_Entropy_Values, EVT_2_MTNR1A_entropy_results$EVT_2_Raw_Entropy_Values, EVT_3_MTNR1A_entropy_results$EVT_3_Raw_Entropy_Values)
# Reshape the combined data into a single column
ALL_EVT_MTNR1A_Raw_entropy <- ALL_EVT_MTNR1A_Raw_entropy %>%
  pivot_longer(everything(), names_to = "EVT_Group", values_to = "EVT_Entropy")

ALL_FB_MTNR1A_Raw_entropy <- bind_cols(FB_1_MTNR1A_entropy_results$FB_1_Raw_Entropy_Values, FB_2_MTNR1A_entropy_results$FB_2_Raw_Entropy_Values, FB_3_MTNR1A_entropy_results$FB_3_Raw_Entropy_Values)
# Reshape the combined data into a single column
ALL_FB_MTNR1A_Raw_entropy <- ALL_FB_MTNR1A_Raw_entropy %>%
  pivot_longer(everything(), names_to = "FB_Group", values_to = "FB_Entropy")

ALL_HB_MTNR1A_Raw_entropy <- bind_cols(HB_1_MTNR1A_entropy_results$HB_1_Raw_Entropy_Values, HB_2_MTNR1A_entropy_results$HB_2_Raw_Entropy_Values, HB_3_MTNR1A_entropy_results$HB_3_Raw_Entropy_Values)
# Reshape the combined data into a single column
ALL_HB_MTNR1A_Raw_entropy <- ALL_HB_MTNR1A_Raw_entropy %>%
  pivot_longer(everything(), names_to = "HB_Group", values_to = "HB_Entropy")

ALL_MIC_MTNR1A_Raw_entropy <- bind_cols(MIC_1_MTNR1A_entropy_results$MIC_1_Raw_Entropy_Values, MIC_2_MTNR1A_entropy_results$MIC_2_Raw_Entropy_Values, MIC_3_MTNR1A_entropy_results$MIC_3_Raw_Entropy_Values)
# Reshape the combined data into a single column
ALL_MIC_MTNR1A_Raw_entropy <- ALL_MIC_MTNR1A_Raw_entropy %>%
  pivot_longer(everything(), names_to = "MIC_Group", values_to = "MIC_Entropy")

ALL_SM_MTNR1A_Raw_entropy <- bind_cols(SM_1_MTNR1A_entropy_results$SM_1_Raw_Entropy_Values, SM_2_MTNR1A_entropy_results$SM_2_Raw_Entropy_Values)
# Reshape the combined data into a single column
ALL_SM_MTNR1A_Raw_entropy <- ALL_SM_MTNR1A_Raw_entropy %>%
  pivot_longer(everything(), names_to = "SM_Group", values_to = "SM_Entropy")

ALL_STB_MTNR1A_Raw_entropy <- bind_cols(STB_1_MTNR1A_entropy_results$STB_1_Raw_Entropy_Values)
# Reshape the combined data into a single column
ALL_STB_MTNR1A_Raw_entropy <- ALL_STB_MTNR1A_Raw_entropy %>%
  pivot_longer(everything(), names_to = "STB_Group", values_to = "STB_Entropy")



####  CALCULATE Z-SCORES ###

ALL_CTB_MTNR1A_entropy_mean <- mean(ALL_CTB_MTNR1A_Raw_entropy$Entropy)
ALL_CTB_MTNR1A_entropy_sd <- sd(ALL_CTB_MTNR1A_Raw_entropy$Entropy)

ALL_CTB_MTNR1A_entropy_Z_scores <- (ALL_CTB_MTNR1A_Raw_entropy$Entropy - ALL_CTB_MTNR1A_entropy_mean)/ALL_CTB_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = ALL_CTB_MTNR1A_entropy_Z_scores)
ALL_CTB_MTNR1A_z_scores<- cbind(ALL_CTB_MTNR1A_Raw_entropy, z_scores_df)
# Overwrite:
write.csv(ALL_CTB_MTNR1A_z_scores, paste0(getwd(), "/Entropy_Results/ALL_CTB_MTNR1A_z_scores.csv"), row.names=FALSE)


ALL_END_MTNR1A_entropy_mean <- mean(ALL_END_MTNR1A_Raw_entropy$Entropy)
ALL_END_MTNR1A_entropy_sd <- sd(ALL_END_MTNR1A_Raw_entropy$Entropy)

ALL_END_MTNR1A_entropy_Z_scores <- (ALL_END_MTNR1A_Raw_entropy$Entropy - ALL_CTB_MTNR1A_entropy_mean)/ALL_CTB_MTNR1A_entropy_sd
# Add z-scores to column in entropy data frame:
z_scores_df <- data.frame(z_scores = ALL_END_MTNR1A_entropy_Z_scores)
ALL_END_MTNR1A_z_scores<- cbind(ALL_END_MTNR1A_Raw_entropy, z_scores_df)
# Overwrite:
write.csv(ALL_END_MTNR1A_z_scores, paste0(getwd(), "/Entropy_Results/ALL_END_MTNR1A_z_scores.csv"), row.names=FALSE)
