# set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary libraries
library(tidyverse)    
library(tidymodels)   
library(limma)      
library(ChAMP)       
library(parallel)   
library(RPMM)
library(glmnet)

# Read training data CSV file and convert it to a data frame
df_train <- read_csv("./Beta_raw_subchallenge1.csv") %>% as.data.frame()
rownames(df_train) <- df_train[,1]    # Set row names to the first column values
df_train[,1] <- NULL                  # Remove the first column

# Read detection P-values CSV file and convert it to a data frame
detp <- read_csv("./DetectionP_subchallenge1.csv") %>% as.data.frame()
rownames(detp) <- detp[,1]            # Set row names to the first column values
detp[,1] <- NULL                      # Remove the first column

# Read sample annotation CSV file
ano <- read.csv("./Sample_annotation.csv")
rownames(ano) <- ano$Sample_ID        # Set row names to the Sample_ID column values
names(ano)[names(ano) == "Sample_ID"] <- "Sample_Name"  # Rename "Sample_ID" column to "Sample_Name"

# Ensure df_train and detp columns match the rownames of ano
df_train <- df_train[, rownames(ano)]
detp <- detp[rownames(df_train), rownames(ano)]

# Filter the beta values using champ.filter
beta_train_filter <- champ.filter(
  beta = df_train,
  pd = ano,
  detP = detp,
  autoimpute = FALSE,
  filterDetP = TRUE,
  filterBeads = FALSE,
  fixOutlier = FALSE,
  filterNoCG = TRUE,
  filterSNPs = TRUE,
  filterMultiHit = TRUE,
  filterXY = TRUE,
  arraytype = "450K"
)$beta

# Store the row names of the filtered beta values
first_filter <- rownames(beta_train_filter)

# Normalize the filtered beta values using BMIQ method
# Define the champ.BMIQ function for normalization of beta values using BMIQ method
champ.BMIQ <- function(beta.v, design.v, nL = 3, doH = TRUE, nfit = 10000, 
                       th1.v = c(0.2, 0.75), th2.v = NULL, niter = 5, tol = 0.001, 
                       sampleID = 1) {
  # Identify indices for type1 and type2 probes based on the design vector
  type1.idx <- which(design.v == 1)
  type2.idx <- which(design.v == 2)
  
  # Extract beta values for type1 and type2 probes
  beta1.v <- beta.v[type1.idx]
  beta2.v <- beta.v[type2.idx]
  
  # Adjust beta values to avoid 0 or 1 by replacing them with the closest non-extreme values
  if (min(beta1.v) == 0) {
    beta1.v[beta1.v == 0] <- min(setdiff(beta1.v, 0))
  }
  if (min(beta2.v) == 0) {
    beta2.v[beta2.v == 0] <- min(setdiff(beta2.v, 0))
  }
  if (max(beta1.v) == 1) {
    beta1.v[beta1.v == 1] <- max(setdiff(beta1.v, 1))
  }
  if (max(beta2.v) == 1) {
    beta2.v[beta2.v == 1] <- max(setdiff(beta2.v, 1))
  }
  
  # Initialize a weight matrix for the EM algorithm
  w0.m <- matrix(0, nrow = length(beta1.v), ncol = nL)
  w0.m[which(beta1.v <= th1.v[1]), 1] <- 1
  w0.m[intersect(which(beta1.v > th1.v[1]), which(beta1.v <= th1.v[2])), 2] <- 1
  w0.m[which(beta1.v > th1.v[2]), 3] <- 1
  
  print("Fitting EM beta mixture to type1 probes")
  set.seed(1234567)  # Set seed for reproducibility
  
  # Randomly sample beta values for fitting the EM model
  rand.idx <- sample(1:length(beta1.v), nfit, replace = FALSE)
  
  # Apply EM algorithm to type1 beta values
  em1.o <- blc(matrix(beta1.v[rand.idx], ncol = 1), w = w0.m[rand.idx, ], maxiter = niter, tol = tol)
  
  # Determine the class of each beta value in type1 probes
  subsetclass1.v <- apply(em1.o$w, 1, which.max)
  
  # Calculate thresholds for class separation in type1 probes
  subsetth1.v <- c(mean(c(max(beta1.v[rand.idx[subsetclass1.v == 1]]), min(beta1.v[rand.idx[subsetclass1.v == 2]]))), 
                   mean(c(max(beta1.v[rand.idx[subsetclass1.v == 2]]), min(beta1.v[rand.idx[subsetclass1.v == 3]]))))
  
  # Assign classes to type1 beta values based on thresholds
  class1.v <- rep(2, length(beta1.v))
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3
  nth1.v <- subsetth1.v
  
  print("Done")
  
  # Compute densities for type1 probe classes
  d1U.o <- density(beta1.v[class1.v == 1])
  d1M.o <- density(beta1.v[class1.v == 3])
  mod1U <- d1U.o$x[which.max(d1U.o$y)]
  mod1M <- d1M.o$x[which.max(d1M.o$y)]
  
  # Compute densities for type2 probes
  d2U.o <- density(beta2.v[which(beta2.v < 0.4)])
  d2M.o <- density(beta2.v[which(beta2.v > 0.6)])
  mod2U <- d2U.o$x[which.max(d2U.o$y)]
  mod2M <- d2M.o$x[which.max(d2M.o$y)]
  
  # Calculate thresholds for class separation in type2 probes
  th2.v <- vector()
  th2.v[1] <- nth1.v[1] + (mod2U - mod1U)
  th2.v[2] <- nth1.v[2] + (mod2M - mod1M)
  
  # Initialize a weight matrix for type2 probes
  w0.m <- matrix(0, nrow = length(beta2.v), ncol = nL)
  w0.m[which(beta2.v <= th2.v[1]), 1] <- 1
  w0.m[intersect(which(beta2.v > th2.v[1]), which(beta2.v <= th2.v[2])), 2] <- 1
  w0.m[which(beta2.v > th2.v[2]), 3] <- 1
  
  print("Fitting EM beta mixture to type2 probes")
  set.seed(1234567)  # Set seed for reproducibility
  
  # Randomly sample beta values for fitting the EM model
  rand.idx <- sample(1:length(beta2.v), nfit, replace = FALSE)
  
  # Apply EM algorithm to type2 beta values
  em2.o <- blc(matrix(beta2.v[rand.idx], ncol = 1), w = w0.m[rand.idx, ], maxiter = niter, tol = tol)
  
  print("Done")
  
  # Determine the class of each beta value in type2 probes
  subsetclass2.v <- apply(em2.o$w, 1, which.max)
  
  # Calculate thresholds for class separation in type2 probes
  subsetth2.v <- c(mean(c(max(beta2.v[rand.idx[subsetclass2.v == 1]]), min(beta2.v[rand.idx[subsetclass2.v == 2]]))), 
                   mean(c(max(beta2.v[rand.idx[subsetclass2.v == 2]]), min(beta2.v[rand.idx[subsetclass2.v == 3]]))))
  
  # Assign classes to type2 beta values based on thresholds
  class2.v <- rep(2, length(beta2.v))
  class2.v[which(beta2.v < subsetth2.v[1])] <- 1
  class2.v[which(beta2.v > subsetth2.v[2])] <- 3
  
  # Calculate average beta values for each class in both type1 and type2 probes
  classAV1.v <- vector()
  classAV2.v <- vector()
  for (l in 1:nL) {
    classAV1.v[l] <- em1.o$mu[l, 1]
    classAV2.v[l] <- em2.o$mu[l, 1]
  }
  
  print("Start normalising type 2 probes")
  
  # Initialize normalized beta values for type2 probes
  nbeta2.v <- beta2.v
  
  # Normalize type2 beta values based on their class
  lt <- 1
  selU.idx <- which(class2.v == lt)
  selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])]
  selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])]
  
  p.v <- pbeta(beta2.v[selUR.idx], em2.o$a[lt, 1], em2.o$b[lt, 1], lower.tail = FALSE)
  q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = FALSE)
  nbeta2.v[selUR.idx] <- q.v
  
  p.v <- pbeta(beta2.v[selUL.idx], em2.o$a[lt, 1], em2.o$b[lt, 1], lower.tail = TRUE)
  q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = TRUE)
  nbeta2.v[selUL.idx] <- q.v
  
  lt <- 3
  selM.idx <- which(class2.v == lt)
  selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])]
  selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])]
  
  p.v <- pbeta(beta2.v[selMR.idx], em2.o$a[lt, 1], em2.o$b[lt, 1], lower.tail = FALSE)
  q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = FALSE)
  nbeta2.v[selMR.idx] <- q.v
  
  # Adjust the values for class 2 if the `doH` parameter is TRUE
  if (doH) {
    lt <- 2
    selH.idx <- c(which(class2.v == lt), selML.idx)
    minH <- min(beta2.v[selH.idx])
    maxH <- max(beta2.v[selH.idx])
    deltaH <- maxH - minH
    deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
    deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])
    nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM
    nminH <- max(nbeta2.v[selU.idx]) + deltaUH
    ndeltaH <- nmaxH - nminH
    hf <- ndeltaH / deltaH
    nbeta2.v[selH.idx] <- nminH + hf * (beta2.v[selH.idx] - minH)
  }
  
  # Combine normalized beta values from type1 and type2 probes
  pnbeta.v <- beta.v
  pnbeta.v[type1.idx] <- beta1.v
  pnbeta.v[type2.idx] <- nbeta2.v
  
  print(paste("Finished for sample ", sampleID, sep = ""))
  
  # Return a list containing the normalized beta values and other relevant information
  return(list(nbeta = pnbeta.v, class1 = class1.v, class2 = class2.v, av1 = classAV1.v, av2 = classAV2.v, hf = hf, th1 = nth1.v, th2 = th2.v))
}

# Load probe information from an RDS file
probeInfoALL.lv <- readRDS("./probe_info_450k.rds")

# Create a design vector based on probe information and the row names of beta_train_filter
design.v <- as.numeric(lapply(probeInfoALL.lv, function(x) x)$Design[match(rownames(beta_train_filter), probeInfoALL.lv$probeID)])

# Apply the champ.BMIQ function to each column of beta_train_filter and normalize the beta values and change the colnames and rownames
beta_norm_BMIQ <- as.data.frame(sapply(1:ncol(beta_train_filter), function(x) champ.BMIQ(beta_train_filter[,x], design.v, sampleID = colnames(beta_train_filter)[x])$nbeta))
rownames(beta_norm_BMIQ) <- rownames(beta_train_filter)
colnames(beta_norm_BMIQ) <- colnames(beta_train_filter)


# Calculate Spearmans' correlations between CpGs and gestational age (GA)
correlations <- apply(as.matrix(beta_norm_BMIQ), 1, function(cpg) cor(as.numeric(cpg), ano[colnames(beta_norm_BMIQ),"GA"], method = "spearman"))

# Combine correlations with normalized beta values, sort by correlations, and remove correlations column
beta_train <- cbind(correlations = correlations, beta_norm_BMIQ)
beta_train <- beta_train[order(beta_train$correlations), ]  
beta_train$correlations <- NULL                             

# Select top and bottom 2000 most correlated CpGs
top_2000 <- head(beta_train, 2000)
bottom_2000 <- tail(beta_train, 2000)
subset_beta_train <- rbind(top_2000, bottom_2000)

# Store the row names of the subset beta values
second_filter <- rownames(subset_beta_train)

# Define a function for feature selection using k-means clustering
feature_selection <- function(subset_beta_train, k){
  
  set.seed(123)
  clustering_result <- kmeans(subset_beta_train, centers = k)  # Perform k-means clustering
  
  subset_beta_train <- cbind(cluster = clustering_result$cluster, subset_beta_train)  # Add cluster labels to data
  
  subset_beta_train_M <- aggregate(. ~ cluster, data = subset_beta_train, FUN = median)  # Aggregate by cluster
  
  cluster_info <- subset_beta_train[, "cluster", drop = FALSE]  # Extract cluster information
  
  return(list(subset_beta_train_M = subset_beta_train_M, cluster_info = cluster_info))
}

# Perform feature selection with k-means clustering
subset_beta_train_M <- feature_selection(subset_beta_train, 1400)
subset_beta_train_mat <- subset_beta_train_M$subset_beta_train_M
subset_beta_train_mat$cluster <- NULL  # Remove cluster column
rownames(subset_beta_train_mat) <- paste("cluster", seq_len(nrow(subset_beta_train_mat)), sep = " ")  # Rename row names

# Extract cluster information
cluster_info <- subset_beta_train_M$cluster_info

# Train a regularized linear regression model using cross-validation
model <- cv.glmnet(x = as.matrix(t(subset_beta_train_mat)), y = log(ano$GA), nfolds = 30, alpha = 0.1, family = "gaussian")

# Save the trained model and other relevant data to files
saveRDS(model, "./model_test_SC1.rds")
saveRDS(first_filter, "./first_filter.rds")
saveRDS(second_filter, "./second_filter.rds")
saveRDS(cluster_info, "./cluster_info.rds")

# Clean up the workspace
rm(list = ls())