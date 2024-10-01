# set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load library
library(optparse)
library(readr)
library(caret)
library(RPMM)
library(glmnet)

# Read the CSV file into a data frame
# The file is located in the directory specified by the --input option
test_data <- as.data.frame(read_csv("./Leaderboard_beta_subchallenge1.csv"))

# Set the row names of the data frame to the values in the first column, and remove the first column from the data frame
rownames(test_data) <- test_data[,1]
test_data <- test_data[,-1]

# Get the Sample IDs for preparing the output file
Sample_IDs <- colnames(test_data)

# Quality control - champ.filter
first_filter <- readRDS("./first_filter.rds")

# Subset the raw data based on the first filter
test_data_subset <- test_data[first_filter, ]

# Normalization for probe 1 and 2 using ChAMP function
# Define the champ.BMIQ function for normalization of beta values using BMIQ method
champ.BMIQ <- function(beta.v, design.v, nL = 3, doH = TRUE, nfit = 10000, 
                       th1.v = c(0.2, 0.75), th2.v = NULL, niter = 5, tol = 0.001, 
                       sampleID = 1) {
  type1.idx <- which(design.v == 1)
  type2.idx <- which(design.v == 2)
  beta1.v <- beta.v[type1.idx]
  beta2.v <- beta.v[type2.idx]
  
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
  
  w0.m <- matrix(0, nrow = length(beta1.v), ncol = nL)
  w0.m[which(beta1.v <= th1.v[1]), 1] <- 1
  w0.m[intersect(which(beta1.v > th1.v[1]), which(beta1.v <= th1.v[2])), 2] <- 1
  w0.m[which(beta1.v > th1.v[2]), 3] <- 1
  
  print("Fitting EM beta mixture to type1 probes")
  set.seed(1234567)
  rand.idx <- sample(1:length(beta1.v), nfit, replace = FALSE)
  em1.o <- blc(matrix(beta1.v[rand.idx], ncol = 1), w = w0.m[rand.idx, ], maxiter = niter, tol = tol)
  
  subsetclass1.v <- apply(em1.o$w, 1, which.max)
  subsetth1.v <- c(mean(c(max(beta1.v[rand.idx[subsetclass1.v == 1]]), min(beta1.v[rand.idx[subsetclass1.v == 2]]))), 
                   mean(c(max(beta1.v[rand.idx[subsetclass1.v == 2]]), min(beta1.v[rand.idx[subsetclass1.v == 3]]))))
  class1.v <- rep(2, length(beta1.v))
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3
  nth1.v <- subsetth1.v
  
  print("Done")
  
  d1U.o <- density(beta1.v[class1.v == 1])
  d1M.o <- density(beta1.v[class1.v == 3])
  mod1U <- d1U.o$x[which.max(d1U.o$y)]
  mod1M <- d1M.o$x[which.max(d1M.o$y)]
  d2U.o <- density(beta2.v[which(beta2.v < 0.4)])
  d2M.o <- density(beta2.v[which(beta2.v > 0.6)])
  mod2U <- d2U.o$x[which.max(d2U.o$y)]
  mod2M <- d2M.o$x[which.max(d2M.o$y)]
  th2.v <- vector()
  th2.v[1] <- nth1.v[1] + (mod2U - mod1U)
  th2.v[2] <- nth1.v[2] + (mod2M - mod1M)
  
  w0.m <- matrix(0, nrow = length(beta2.v), ncol = nL)
  w0.m[which(beta2.v <= th2.v[1]), 1] <- 1
  w0.m[intersect(which(beta2.v > th2.v[1]), which(beta2.v <= th2.v[2])), 2] <- 1
  w0.m[which(beta2.v > th2.v[2]), 3] <- 1
  
  print("Fitting EM beta mixture to type2 probes")
  set.seed(1234567)
  rand.idx <- sample(1:length(beta2.v), nfit, replace = FALSE)
  em2.o <- blc(matrix(beta2.v[rand.idx], ncol = 1), w = w0.m[rand.idx, ], maxiter = niter, tol = tol)
  
  print("Done")
  
  subsetclass2.v <- apply(em2.o$w, 1, which.max)
  subsetth2.v <- c(mean(c(max(beta2.v[rand.idx[subsetclass2.v == 1]]), min(beta2.v[rand.idx[subsetclass2.v == 2]]))), 
                   mean(c(max(beta2.v[rand.idx[subsetclass2.v == 2]]), min(beta2.v[rand.idx[subsetclass2.v == 3]]))))
  class2.v <- rep(2, length(beta2.v))
  class2.v[which(beta2.v < subsetth2.v[1])] <- 1
  class2.v[which(beta2.v > subsetth2.v[2])] <- 3
  
  classAV1.v <- vector()
  classAV2.v <- vector()
  for (l in 1:nL) {
    classAV1.v[l] <- em1.o$mu[l, 1]
    classAV2.v[l] <- em2.o$mu[l, 1]
  }
  
  print("Start normalising type 2 probes")
  nbeta2.v <- beta2.v
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
  
  pnbeta.v <- beta.v
  pnbeta.v[type1.idx] <- beta1.v
  pnbeta.v[type2.idx] <- nbeta2.v
  
  print(paste("Finished for sample ", sampleID, sep = ""))
  return(list(nbeta = pnbeta.v, class1 = class1.v, class2 = class2.v, av1 = classAV1.v, av2 = classAV2.v, hf = hf, th1 = nth1.v, th2 = th2.v))
}

# Load probe information from an RDS file
probeInfoALL.lv <- readRDS("./probe_info_450k.rds")

# Create a design vector based on probe information and the row names of beta_train_filter
design.v <- as.numeric(lapply(probeInfoALL.lv, function(x) x)$Design[match(rownames(test_data_subset), probeInfoALL.lv$probeID)])

# Apply the champ.BMIQ function to each column of beta_train_filter and normalize the beta values and change the colnames and rownames
test_data_subset_norm <- as.data.frame(sapply(1:ncol(test_data_subset), function(x) champ.BMIQ(test_data_subset[,x], design.v, sampleID = colnames(test_data_subset)[x])$nbeta))
rownames(test_data_subset_norm) <- rownames(test_data_subset)
colnames(test_data_subset_norm) <- colnames(test_data_subset)

# Correlation filter
second_filter <- readRDS("./second_filter.rds")

# Subset the normalized data based on the second filter
test_data_subset_norm <- test_data_subset_norm[second_filter,]

# Clustering information to cluster the filtered normalized data
cluster_info <- readRDS("./cluster_info.rds")

# Function for feature selection/engineering using k-means clustering
# Steps : Perform k-means clustering, Combine cluster information with subset data, Calculate median values for each cluster, Extract cluster information
# Rename the rows of the selected features

test_data_subset_norm <- test_data_subset_norm[match(rownames(cluster_info), rownames(test_data_subset_norm)), ]
test_data_subset_norm <- cbind(cluster = cluster_info, test_data_subset_norm)
test_data_subset_norm <- aggregate(. ~ cluster, data = test_data_subset_norm, FUN = median)
test_data_subset_norm$cluster <- NULL
rownames(test_data_subset_norm) <- paste("cluster", seq_len(nrow(test_data_subset_norm)), sep = " ")

# Load the pre-trained model from an RDS file
model <- readRDS("./model_test_SC1.rds")

# Transpose the selected feature matrix
test_data <- as.data.frame(t(test_data_subset_norm))

# Make predictions using the pre-trained model
ga <- as.numeric(predict(model, newx = as.matrix(test_data), s = "lambda.min"))
ga <- exp(ga)

# Ensure predictions are within the valid range [5, 44]
ga[ga > 44] <- 44
ga[ga < 5] <- 5

# Prepare the output data frame with sample IDs and predictions
output_df <- data.frame(
  Sample = Sample_IDs,
  GA_prediction = ga
)

# Rename the columns of the output data frame
names(output_df) <- c("ID", "GA_prediction")

# Write the output data frame to a CSV file in the specified output directory
write_csv(output_df, "./predictions.csv")

