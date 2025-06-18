######################################
#### This code runs the simulation study in
#### High-Dimensional Block Diagonal Covariance Structure Detection Using Singular Vectors
#### published in the Journal of Computational and Graphical Statistics

library(bdsvd)    #for BD-SVD
library(cvCovEst) #for the ad hoc procedure
library(mvtnorm)  #for rmvnorm
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

design = "c" #Choose simulation design "a", "b", "c", or "d"

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

for (s in 1:S){
  set.seed(123+s)
  
  #Simulate data matrix X
  Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
  X <- rmvnorm(n, mean=rep(0, p), sigma=Sigma)
  colnames(X) = 1:ncol(X)
  
  
  
  
  if(BDSVD){
    cat("Compute BD-SVD\r")
    
    start <- proc.time()
    SVBD <- bdsvd(X, standardize = FALSE, trace = TRUE)
    SVBD.TIME[s] <- (proc.time() - start)[3]/60
    
    ConfusionMatrix <- get.confusionMatrix(blocks = SVBD, labels = labels, b = b)
    SVBD.SENS[s] <- ConfusionMatrix$sensitivity
    SVBD.SPEC[s] <- ConfusionMatrix$specificity
    SVBD.FDR[s]  <- ConfusionMatrix$FDR
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


Results <- cbind.data.frame(SVBD.SENS, SVBD.SPEC, SVBD.FDR, SVBD.TIME,
                            Est.1.SENS, Est.1.SPEC, Est.1.FDR, Est.1.TIME,
                            Est.2.SENS, Est.2.SPEC, Est.2.FDR, Est.2.TIME,
                            SHDJ.SENS, SHDJ.SPEC, SHDJ.FDR, SHDJ.TIME,
                            SHRR.SENS, SHRR.SPEC, SHRR.FDR, SHDJ.TIME
)
head(Results)


# ─────────────────────────────────────────────────────────────────
# ► SAVE & DISPLAY SYSTEM FOR ALL 3 METHODS
# ─────────────────────────────────────────────────────────────────

library(tidyr)
library(ggplot2)

# 1) Build long data.frames for each metric across BD-SVD, Est.1, and Est.2
df_sens <- data.frame(
  `BD-SVD` = SVBD.SENS,
  Est.1    = Est.1.SENS,
  Est.2    = Est.2.SENS
) %>% pivot_longer(everything(), names_to="Method", values_to="Sensitivity")

df_spec <- data.frame(
  `BD-SVD` = SVBD.SPEC,
  Est.1    = Est.1.SPEC,
  Est.2    = Est.2.SPEC
) %>% pivot_longer(everything(), names_to="Method", values_to="Specificity")

df_fdr <- data.frame(
  `BD-SVD` = SVBD.FDR,
  Est.1    = Est.1.FDR,
  Est.2    = Est.2.FDR
) %>% pivot_longer(everything(), names_to="Method", values_to="FDR")


# 2) Create the three boxplots
p_sens <- ggplot(df_sens, aes(x = Method, y = Sensitivity)) +
  geom_boxplot() +
  ylim(0, 1) +
  ggtitle(sprintf("Sensitivity (n=%d, p=%d, b=%d, design='%s')", n, p, b, design)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )

p_spec <- ggplot(df_spec, aes(x = Method, y = Specificity)) +
  geom_boxplot() +
  ylim(0, 1) +
  ggtitle(sprintf("Specificity (n=%d, p=%d, b=%d, design='%s')", n, p, b, design)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )

p_fdr <- ggplot(df_fdr, aes(x = Method, y = FDR)) +
  geom_boxplot() +
  ylim(0, 1) +
  ggtitle(sprintf("FDR (n=%d, p=%d, b=%d, design='%s')", n, p, b, design)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )


# 3) Dynamic filename base
base_fn <- sprintf("n%d_p%d_b%d_design%s", n, p, b, design)

# 4) Save each plot as its own PNG
ggsave(filename = paste0(base_fn, "_sensitivity.png"),
       plot     = p_sens, width = 5, height = 5, dpi = 300)

ggsave(filename = paste0(base_fn, "_specificity.png"),
       plot     = p_spec, width = 5, height = 5, dpi = 300)

ggsave(filename = paste0(base_fn, "_fdr.png"),
       plot     = p_fdr,  width = 5, height = 5, dpi = 300)

message("Saved plots: ",
        paste0(base_fn, c("_sensitivity.png",
                          "_specificity.png",
                          "_fdr.png"), collapse = ", "))




# 6) Save entire workspace for later
save.image(file = paste0(base_fn, ".Rdata"))
message("Workspace saved to ", base_fn, ".Rdata")






# 5) Display them in R’s plotting window
print(p_sens)
print(p_spec)
print(p_fdr)

