######################################
#### This code runs the simulation study in
#### High-Dimensional Block Diagonal Covariance Structure Detection Using Singular Vectors
#### published in the Journal of Computational and Graphical Statistics

# Source path to BDSVD-FINAL.R
source("PATH/GitHub/Thesis/BDSVD-FINAL.R")
library(cvCovEst) #for the ad hoc procedure
library(mvtnorm)  #for rmvnorm1
library(shock)    #for SHDJ and SHRR

get.confusionMatrix <- function(blocks, labels, b){
  TP <- 0
  FN <- 0
  FP <- 0
  TN <- 0
  
  block.labels <- list()
  for(i in 1:length(blocks)){
    block.labels[[i]] <- as.integer(blocks[[i]])
  }
  
  SigmaEst <- matrix(FALSE, p, p)
  for(i in 1:length(blocks)){
    SigmaEst[block.labels[[i]], block.labels[[i]] ] <- TRUE
  }
  
  for(i in 1:b){
    if(all(SigmaEst[which(labels == i), which(labels == i)])){
      TN <- TN + 1 #Not splitting within Sigma_i
    } else{
      FP <- FP + 1 #Splitting within Sigma_i
    }
  }
  
  for(i in 1:b){
    if(any(SigmaEst[which(labels != i), which(labels == i)])){
      FN <- FN + 1 #Not splitting Sigma_i and Sigma_j
    } else{
      TP <- TP + 1 #Spliting Sigma_i and Sigma_j
    }
  }
  
  sens <- TP / (TP + FN)
  spec <- TN / (FP + TN)
  FPR <- FP / (FP + TN)
  FDR <- FP / (TP + FP)
  
  if(TP+FP == 0){FDR <- 0}
  
  return(list(TP=TP, TN=TN, FP=FP,FN=FN,sensitivity=sens,specificity=spec,FDR=FDR, FPR=FPR))
}

p <- 500 #number of variables
n <- 250 #number of observations
b <- 10  #number of blocks
S <- 100 #Number of simulations

design = "c" #Choose simulation design "a", "b", "c", "d" or "e"


# --- BD-SVD Algorithm Parameters ---
# These are used inside the loop to configure each run of the `bdsvd` function.
q.selection <- TRUE  # TRUE: Automatically select 'q' at each split. FALSE: Use q=1.
mode <- "ERM"        # Method for q-selection. Options: "ERM", "EKC", "PAR".
weight_method <- "equal" # Sparsity allocation for q>1 splits. Options: "equal", "explained_var".



BDSVD <- TRUE #Perform BD-SVD? TRUE or FALSE
AdHoc <- TRUE #Perform the ad hoc procedure? TRUE or FALSE
SH    <- FALSE #Perform SHDJ and SHRR? TRUE or FALSE

labels <- factor(rep(1:b, length.out=p))
labels <- labels[order(labels)]

SVBD.SENS <- vector(length=S)
SVBD.SPEC <- vector(length=S)
SVBD.FDR  <- vector(length=S)
SVBD.TIME <- vector(length=S)

SHDJ.SENS <- vector(length=S)
SHDJ.SPEC <- vector(length=S)
SHDJ.FDR  <- vector(length=S)
SHDJ.TIME <- vector(length=S)

SHRR.SENS <- vector(length=S)
SHRR.SPEC <- vector(length=S)
SHRR.FDR  <- vector(length=S)

Est.1.SENS <- vector(length=S)
Est.1.SPEC <- vector(length=S)
Est.1.FDR  <- vector(length=S)
Est.1.TIME <- vector(length=S)

Est.2.SENS <- vector(length=S)
Est.2.SPEC <- vector(length=S)
Est.2.FDR  <- vector(length=S)
Est.2.TIME <- vector(length=S)


q_count <- integer(S)
all_qs <- integer(0)


for (s in 1:S){
  set.seed(123+s)
  
  #Simulate data matrix X
  Sigma <- bdsvd.cov.sim(p = p, b = b, design = design, noise_level = 0.15) #set noise level for design "e"
  X <- rmvnorm(n, mean=rep(0, p), sigma=Sigma)
  colnames(X) = 1:ncol(X)
  
  if(BDSVD){
    cat("Compute BD-SVD\r")
    
    start <- proc.time()
    SVBD <- bdsvd(X , q.selection = q.selection, mode = mode, standardize = FALSE,
                  trace = TRUE, weight_method = weight_method)
    
    SVBD.TIME[s] <- (proc.time() - start)[3]/60
    
    ConfusionMatrix <- get.confusionMatrix(blocks = SVBD, labels = labels, b = b)
    SVBD.SENS[s] <- ConfusionMatrix$sensitivity
    SVBD.SPEC[s] <- ConfusionMatrix$specificity
    SVBD.FDR[s]  <- ConfusionMatrix$FDR
    
    
    q_count[s] <- sum(attr(SVBD, "chosen_qs") > 1)
    all_qs <- c(all_qs, attr(SVBD, "chosen_qs"))
  }
  
  
  
  
  if(AdHoc){
    cat("Compute Est.1 and Est.2\r")
    
    #Est.1
    start <- proc.time()
    blocks.1 <- detect.blocks(scadEst(X, 0.2), 0.1)
    Est.1.TIME[s] <- (proc.time() - start)[3]/60
    
    Estblocks <- list()
    for(i in 1:length(blocks.1)){
      Estblocks[[i]] <- as.character(blocks.1[[i]]@features)
    }
    
    ConfusionMatrix <- get.confusionMatrix(blocks = Estblocks, labels = labels, b = b)
    Est.1.SENS[s] <- ConfusionMatrix$sensitivity
    Est.1.SPEC[s] <- ConfusionMatrix$specificity
    Est.1.FDR[s]  <- ConfusionMatrix$FDR
    
    #Est.2
    start <- proc.time()
    blocks.2 <- detect.blocks(scadEst(X, 0.2), 0.2)
    Est.2.TIME[s] <- (proc.time() - start)[3]/60
    
    Estblocks <- list()
    for(i in 1:length(blocks.2)){
      Estblocks[[i]] <- as.character(blocks.2[[i]]@features)
    }
    
    ConfusionMatrix <- get.confusionMatrix(blocks = Estblocks, labels = labels, b = b)
    Est.2.SENS[s] <- ConfusionMatrix$sensitivity
    Est.2.SPEC[s] <- ConfusionMatrix$specificity
    Est.2.FDR[s]  <- ConfusionMatrix$FDR
  }
  
  
  
  
  if(SH){
    cat("Compute SHDJ and SHRR\r")
    
    start <- proc.time()
    resShock <- shockSelect(X)
    SHDJ.TIME[s] <- (proc.time() - start)[3]/60
    
    ## SHDJ
    shdjAEstim <- diag(p)
    for(i in 1:length(unique(resShock$SHDJlabels))){
      
      stepdata <- as.matrix(X[,resShock$SHDJlabels==i],nrow=dim(X)[1])
      if(dim(stepdata)[2]>1){
        resNet <- networkInferenceGlassoBIC(stepdata)
        shdjAEstim[resShock$SHDJlabels==i,resShock$SHDJlabels==i] <- resNet$A
      }
    }
    shdjAEstim.blocks <- list()
    for(i in 1:length(detect.blocks(shdjAEstim, 0))){
      shdjAEstim.blocks[[i]] <- as.character(detect.blocks(shdjAEstim, 0)[[i]]@features)
    }
    
    ConfusionMatrix <- get.confusionMatrix(blocks = shdjAEstim.blocks, labels = labels, b = b)
    SHDJ.SENS[s] <- ConfusionMatrix$sensitivity
    SHDJ.SPEC[s] <- ConfusionMatrix$specificity
    SHDJ.FDR[s]  <- ConfusionMatrix$FDR
    
    #SHRR
    shrrAEstim <- diag(p)
    for(i in 1:length(unique(resShock$SHRRlabels))){
      
      stepdata <- as.matrix(X[,resShock$SHRRlabels==i],nrow=dim(X)[1])
      if(dim(stepdata)[2]>1){
        resNet <- networkInferenceGlassoBIC(stepdata)
        shrrAEstim[resShock$SHRRlabels==i,resShock$SHRRlabels==i] <- resNet$A
      }
    }
    shrrAEstim.blocks <- list()
    for(i in 1:length(detect.blocks(shrrAEstim, 0))){
      shrrAEstim.blocks[[i]] <- as.character(detect.blocks(shrrAEstim, 0)[[i]]@features)
    }
    
    ConfusionMatrix <- get.confusionMatrix(blocks = shrrAEstim.blocks, labels = labels, b = b)
    SHRR.SENS[s] <- ConfusionMatrix$sensitivity
    SHRR.SPEC[s] <- ConfusionMatrix$specificity
    SHRR.FDR[s]  <- ConfusionMatrix$FDR
  }
  
  
  print(paste0("======================="))
  print(paste0("== RESULT FOR s in 1:", s))
  print(paste0("======================="))
  print(paste0("AVERAGE SENSITIVITY:"))
  print(paste0("BD-SVD:",  round(mean(SVBD.SENS[1:s]), 2)  ))
  print(paste0("Est.1: ",  round(mean(Est.1.SENS[1:s]), 2) ))
  print(paste0("Est.2: ",  round(mean(Est.2.SENS[1:s]), 2) ))
  print(paste0("SHDJ:  ",  round(mean(SHDJ.SENS[1:s]), 2)  ))
  print(paste0("SHRR:  ",  round(mean(SHRR.SENS[1:s]), 2)  ))
  print(paste0("======================="))
  print(paste0("AVERAGE SPECIFICITY:"))
  print(paste0("BD-SVD:",  round(mean(SVBD.SPEC[1:s]), 2)  ))
  print(paste0("Est.1: ",  round(mean(Est.1.SPEC[1:s]), 2) ))
  print(paste0("Est.2: ",  round(mean(Est.2.SPEC[1:s]), 2) ))
  print(paste0("SHDJ:  ",  round(mean(SHDJ.SPEC[1:s]), 2)  ))
  print(paste0("SHRR:  ",  round(mean(SHRR.SPEC[1:s]), 2)  ))
  print(paste0("======================="))
  print(paste0("AVERAGE FDR:"))
  print(paste0("BD-SVD:",  round(mean(SVBD.FDR[1:s]), 2)  ))
  print(paste0("Est.1: ",  round(mean(Est.1.FDR[1:s]), 2) ))
  print(paste0("Est.2: ",  round(mean(Est.2.FDR[1:s]), 2) ))
  print(paste0("SHDJ:  ",  round(mean(SHDJ.FDR[1:s]), 2)  ))
  print(paste0("SHRR:  ",  round(mean(SHRR.FDR[1:s]), 2)  ))
  print(paste0("======================="))
  print(paste0("AVERAGE TIME (min.):"))
  print(paste0("BD-SVD:",  round(mean(SVBD.TIME[1:s]), 2)  ))
  print(paste0("Est.1: ",  round(mean(Est.1.TIME[1:s]), 2) ))
  print(paste0("Est.2: ",  round(mean(Est.2.TIME[1:s]), 2) ))
  print(paste0("SHDJ:  ",  round(mean(SHDJ.TIME[1:s]), 2)  ))
  print(paste0("SHRR:  ",  round(mean(SHDJ.TIME[1:s]), 2)  ))
  print(paste0("======================="))
}


Results <- cbind.data.frame(SVBD.SENS, SVBD.SPEC, SVBD.FDR, SVBD.TIME,q_count,
                            Est.1.SENS, Est.1.SPEC, Est.1.FDR, Est.1.TIME,
                            Est.2.SENS, Est.2.SPEC, Est.2.FDR, Est.2.TIME,
                            SHDJ.SENS, SHDJ.SPEC, SHDJ.FDR, SHDJ.TIME,
                            SHRR.SENS, SHRR.SPEC, SHRR.FDR, SHDJ.TIME
)
head(Results)

library(tidyr); library(dplyr); library(ggplot2); library(patchwork)

# --- Results Processing, Plotting, and Saving ---

tryCatch({
  
  # 1. Create a data frame ONLY for the BD-SVD results
  df_bdsvd <- data.frame(
    Sensitivity = SVBD.SENS,
    Specificity = SVBD.SPEC,
    FDR = SVBD.FDR
  ) %>%
    tidyr::pivot_longer(everything(), names_to = "Metric", values_to = "Value")
  
  # 2. Create the individual boxplots (now black and white)
  p_sens <- ggplot(df_bdsvd %>% filter(Metric == "Sensitivity"),
                   aes(x = Metric, y = Value)) +
    geom_boxplot(fill = "white", color = "black") +
    ylim(0, 1) +
    ggtitle("Sensitivity") +
    theme_bw() +                                      # ← white panel + grey grid
    theme(
      panel.grid.minor = element_blank(),             # (optional) drop minor grid
      plot.title       = element_text(hjust = 0.5),
      axis.title.x     = element_blank(),
      axis.text.x      = element_blank()
    )
  
  p_spec <- ggplot(df_bdsvd %>% filter(Metric == "Specificity"),
                   aes(x = Metric, y = Value)) +
    geom_boxplot(fill = "white", color = "black") +
    ylim(0, 1) +
    ggtitle("Specificity") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(hjust = 0.5),
      axis.title.x     = element_blank(),
      axis.text.x      = element_blank()
    )
  
  p_fdr <- ggplot(df_bdsvd %>% filter(Metric == "FDR"),
                  aes(x = Metric, y = Value)) +
    geom_boxplot(fill = "white", color = "black") +
    ylim(0, 1) +
    ggtitle("FDR") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(hjust = 0.5),
      axis.title.x     = element_blank(),
      axis.text.x      = element_blank()
    )
  
  # 3. Define dynamic titles and filenames
  title_text   <- sprintf(
    "BD-SVD Results for Design '%s' (n=%d, p=%d, b=%d)",
    design, n, p, b
  )
  caption_text <- sprintf(
    "Mode: '%s' | Weight Method: '%s'",
    mode, weight_method
  )
  base_filename <- sprintf(
    "n%d_p%d_b%d_design%s_mode%s_weight%s",
    n, p, b, design, mode, weight_method
  )
  
  plot_filename <- paste0(base_filename, ".png")
  q_hist_filename <- paste0(base_filename, "_q_histogram.png")
  
  # 4. Combine plots into a single figure
  final_plot <- (p_sens | p_spec | p_fdr) +
    plot_annotation(
      title = title_text,
      caption = caption_text
    ) &
    theme(plot.caption = element_text(hjust = 0, face = "italic"))
  
  # 5. Save the plot
  ggsave(plot_filename, plot = final_plot, width = 8, height = 4, dpi = 300)
  print(paste("Plot saved as:", plot_filename))
  
  # 6. Save the q-histogram (now black and white)
  if (length(all_qs) > 0) {
    png(q_hist_filename)
    hist(all_qs,
         breaks = seq(min(all_qs, na.rm = TRUE) - 0.5, max(all_qs, na.rm = TRUE) + 0.5, by = 1),
         main   = "Histogram of Chosen q Values",
         xlab   = "q",
         ylab   = "Frequency",
         col    = "white",  # MODIFIED: for B&W
         border = "black") # MODIFIED: for B&W
    dev.off()
    print(paste("Histogram of q values saved as:", q_hist_filename))
  } else {
    print("No q > 1 splits occurred, so no q-histogram was generated.")
  }
  
}, error = function(e) {
  cat("An error occurred during plotting/saving figures: ", e$message, "\n")
  cat("Skipping plot generation and proceeding to save data.\n")
})

# 7. Save the raw results to an .Rdata file
data_filename <- sprintf("n%d_p%d_design%s_mode%s_weight%s.Rdata",
                         n, p, design, mode, weight_method)

save.image(file = data_filename)

print(paste("Complete workspace saved to:", data_filename))

