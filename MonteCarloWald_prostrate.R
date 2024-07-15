####Monte Carlo Wald test for Prostrate Cancer dataset
### Similar to the whole EM algotihm but recompute sampling probabilities and weights at 
### May 5th Update:
library(dplyr)

#----#---#----#----# The EM algorithm #----#----#----#-----#-----#-----#
## Modify the function to take in two datafiles: Complete dataset and Kinship matrix:
## EPS data will be generated in the function:
## Read in thetahat file and other files
## Do this for each SNP. Insert standard error computation for beta and betag
## Edit the 
## EPS for real data with some level of missing?? How to deal with the missing genotypes? 
## locate it and add to the EPS portion. 

##Program takes in data, kinship matrix and already computed theta hat for each SNP. 


Em_Wald_test<- function(data, M=100,  K, theta1){
  
  ##Define missing genotypes
  geno_comp<- data[,3]
  Xstar<- as.matrix(cbind(x0=1, x1=geno_comp)) # Complete genotype vector with intercept
  
  #First obtain the EPS data: 
  ## For simulated data: we want to retain the top and bottom 20%
  n_rows<- nrow(data) ##already sorted
  ##Top and bottom indices
  top_20_percent_index <- round(0.2 * n_rows)
  bottom_20_percent_index <- (n_rows - round(0.2 * n_rows)) + 1
  
  top_20_percent<- data[1:top_20_percent_index,]
  bottom_20_percent<- data[bottom_20_percent_index:n_rows,]
  #k<- rbind(top_20_percent, bottom_20_percent) ##or
  to_keep<- c(1:top_20_percent_index, bottom_20_percent_index:n_rows)
  EPS_data<- data[to_keep, ]
  ##Now the missing data are the ones obtained by substracting (data - EPS_data)
  dat_missing<- subset(data, !(data$ID %in%  EPS_data$ID))
  
  ## Observed genotype frequencies to be used for initial starting point for the algorithm
  genotype_counts<- as.data.frame(table(geno_comp))  
  geno_prop<- c()
  for (i in genotype_counts["Freq"]){
    d<-length(geno_comp)
    geno_prop<- i/d
  }
  ### Read in theta hat file 
  # theta0<- c()
  # theta0[1:ncol(Xstar)]<- c(-0.0523, 0.0477)
  # theta0[3]<- 0.3  ###sigmasqub
  # theta0[4]<- 1.001904 ##sigmasque
  # theta0[5:7] <- geno_prop
  
  theta1<- c()
  theta1[1:ncol(Xstar)]<- c(-0.0523, 0.0477)
  theta0[3]<- 0.3  ###sigmasqub
  theta0[4]<- 1.001904 ##sigmasque
  theta0[5:7] <- geno_prop
  
  id <- data[,1]
  u.id <- unique(id)
  n <- length(u.id)
  y<- data[,2]
  #K<- kin_mat
  beta0<- theta0[1:ncol(Xstar)]
  sigmasqub0<- theta0[3] ##Variance component for the randomeffect
  sigmasque0<- theta0[4] ## Variance component of the residual error
  initial.theta <- theta0
  theta.hist <- theta0
  
  #for (it in 1:niter) {
  D <- (sigmasqub0 * K)
  invD<- Matrix::chol2inv(D)
  ## Variance of y
  idenRand<- diag(sigmasque0, nrow=n, ncol=n)
  Vary<- D + idenRand
  ## variance due to the suffcient statistics: check again
  Vy<- Matrix::chol2inv(invD + diag(rep(sigmasque0, n))) 
  #Vy<- solve((invD + sigmasque0), tol = 1e-17)
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----
  #STEP2: Compute weights for the individuals with missing genotype data: 
  ## These are the individuals in the middle of the distribution after 
  ## we have obtained the EPS samples: The data set is dat_missing
  #-----||-----||-----||-----||-----||-----||-----||-----||----
  #-----||-----||-----||-----||-----||-----||-----||-----||----
  
  y_mis<- dat_missing$BMI
  id_mis<-dat_missing$BMI
  missing_geno_weights<- matrix(rep(0, length(y_mis) * 3), ncol = 3)
  Ind_var<- Vary[row(Vary) == col(Vary)] ###diagonal elements of the variance matrix signifying   
  Missing_ind_var<- Ind_var[as.integer(rownames(dat_missing))]
  
  ## Weights for the missing individuals
  miss_geno_probs<- matrix(rep(0, length(y_mis) * 3), ncol = 3)
  for (i in seq_along(y_mis)) {
    for (j in 0:2){
      #y_imis <- y_mis[id_mis == u.idMis[i]]
      miss_geno_probs[i, j+1]<- dnorm(y_mis[i], (theta0[1] + j*theta0[2]), sqrt(Missing_ind_var[i])) * geno_prop[j+1]
    }
  }
  
  missing_geno_weights<- matrix(rep(0, length(y_mis) * 3), ncol = 3)
  for (i in seq_along(y_mis)){
    for (j in seq_along(1:3)){
      missing_geno_weights[i,j]<- miss_geno_probs[i,j]/sum(miss_geno_probs[i,])
    }
  }
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----
  #STEP2b: Prepare the weights for the individuals with complete genotype data 
  # Recall the weight is 1 for the actual genotype of the individual and 0 for the
  # other two genotypes. We show this below:
  #-----||-----||-----||-----||-----||-----||-----||-----||----
  #-----||-----||-----||-----||-----||-----||-----||-----||----
  
  Weights_complete<- matrix(rep(0, nrow(EPS_data) * 3), ncol = 3)
  for (i in 1:nrow(EPS_data)){
    if (EPS_data$Geno[i] =='0'){
      Weights_complete[i,1]<- 1
      Weights_complete[i,2]<- 0
      Weights_complete[i,3]<- 0
    } else if (EPS_data$Geno[i] =='1'){
      Weights_complete[i,1]<- 0
      Weights_complete[i,2]<- 1
      Weights_complete[i,3]<- 0
    }
    else{
      Weights_complete[i,1]<- 0
      Weights_complete[i,2]<- 0
      Weights_complete[i,3]<- 1
    }
  }
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----
  #STEP2c: Monte Carlo Sampling:
  # A genotype with level 0,1,2 is sampled for each individual 
  # with missing genotype based vectors from the distribution.
  # missing_geno_weights.
  #-----||-----||-----||-----||-----||-----||-----||-----||----
  #-----||-----||-----||-----||-----||-----||-----||-----||----
  M=100
  sample_outcomes<- matrix(ncol=nrow(missing_geno_weights), nrow = M)
  genof= c(0,1,2)
  for (mis_ind in 1:nrow(missing_geno_weights)){
    for (l in 1:M){
      sample_outcomes[l,mis_ind]<- sample(genof, 1, replace = TRUE, prob = missing_geno_weights[mis_ind,])
    }
  }
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2d: Obtain the g_l i.e the lth genotype vector over all individuals.
  # This is a combination of gl= g_obs, g_(mis,l) where g_obs is the vector of 
  # complete genotypes and g_(mis_l) is the lth sampled genotype vector for the 
  # missing individuals. I want to extract every line of the matrix 
  # sample_outcomes and combine it with g_obs to form g_l.
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  
  g_obs<- c(EPS_data$Geno)
  gl <-list()
  for (i in 1:nrow(sample_outcomes)){
    gl[[i]] <- c(g_obs, sample_outcomes[i,])
  }
  
  #-----//----------//----------//--------
  ## Now for each lth sampled genotype vector, combine it with the intercept column
  X_lsampled <- list()
  for (i in 1:nrow(sample_outcomes)){
    X_lsampled[[i]]<-   as.matrix(cbind(Xstar[,1], gl[[i]]))
  }
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2e: Obtain b_l and V_l at the t+1th iteration
  # At the t=0 iteration: 
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  
  gll  <- list()
  b_l1 <- list()
  for (l in 1:M){
    gll[[l]]<- (y - (X_lsampled[[l]]%*%matrix(beta0)))
    b_l1[[l]] <-  (1/sigmasque0) * Vy %*% gll[[l]]
  }
  
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2f:Compute betastar at the t+1th iteration.
  # At the t=0 iteration: 
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  
  D<- list()
  for (l in 1:M){
    D[[l]]<- t(X_lsampled[[l]]) %*% X_lsampled[[l]]
  }
  XX<-solve(Reduce('+', D))
  
  F1<- list()
  F2<- list()
  for (l in 1:M){
    F1[[l]]<- t(y - b_l1[[l]])
    F2[[l]] <- F1[[l]] %*% X_lsampled[[l]]
  }
  XZ<- t(Reduce('+', F2))
  betastar1<- (XX %*% XZ)
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2g:Estimating covariance parameters at the t+1th iteration.
  # At the t=0 iteration: 
  # For sigmaque
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  #E<-0
  ee <- c()
  e<- list()
  #ee<- list()
  for (l in 1:M){
    e[[l]] <- (y - (X_lsampled[[l]] %*% betastar1) - b_l1[[l]])
    #ee[[l]]<- t(e[[l]]) %*% e[[l]]
    ee<- sum(t(e[[l]]) %*% e[[l]])
    
  }
  # E<- sum(ee)
  #E<- Reduce('+', ee)
  
  #sigmasque1<-(1/n) * (sum(diag(Vy))) + ((1/M) * ee)
  sigmasque1<- (1/n) * ( (sum(diag(Vy))) + ((1/M) * ee))
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2f:Estimating covariance parameters at the t+1th iteration.
  # At the t=0 iteration: 
  # For sigmasqub
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  #B<-0
  invK<- solve(K)
  bk<- list()
  for (l in seq_along(M)){
    bk[[l]]<- t(b_l1[[l]]) %*% invK %*% b_l1[[l]]
  }
  B <- Reduce('+', bk)
  kv<- sum(diag(solve(K) %*% Vy)) ##trace of k inverse and vy
  sigmasqub1<- (1/n) * ( kv + (1/M)*(B))
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2g:Estimating genetic covariate parameters at the t+1th iteration.
  # At the t=0 iteration: 
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  Nums<- nrow(Weights_complete)
  to_keep<- c(1:top_20_percent_index, bottom_20_percent_index:n_rows)
  w1<- Weights_complete[c(1:top_20_percent_index), ]
  w3<- Weights_complete[c((top_20_percent_index + 1):Nums),]
  w2<- missing_geno_weights
  All_ind_weights<- as.data.frame(rbind(w1,w2,w3))
  Obs_freq1<- length(which(All_ind_weights$V1 =="1"))
  miss_freq1<- All_ind_weights$V1[All_ind_weights$V1 < 1]
  gamma0<-  (Obs_freq1/n) + (sum(miss_freq1)/n)
  gamma1<- (length(which(All_ind_weights$V2 =="1"))/n) + sum(All_ind_weights$V2[All_ind_weights$V2 <1])/n
  gamma2<- 1-gamma0-gamma1
  geno_prop1 <- c(gamma0,gamma1,gamma2)
  
  
  thetahat <- c(betastar1, sigmasqub1, sigmasque1,geno_prop1)
  #theta.hist <- rbind(theta.hist, theta1)
  
  
  ### Monte Carlo approach to the Wald test.
  ### Second derivative for betag
  ## Gives a two by two matrix: To resuse this in the information matrix,
  ## we pick the diagonal elements which correspond to the second derivatives of
  ## the intercept and betag
  ##Convert the theta1 to vector format for use in subsequent computations
  t1<- c(theta1[,1], theta1[,2], theta1[,3], theta1[,4])
  t2 <- t1[1:2]
  
  Xll<- Reduce('+', D) ##
  betagsqu<- (1/(M * t1[3])) * Xll ##compute this after obtaining thetahat
  
  ### First derivative computation:
  ##Initialize list to store components of the bracket
  FF1<- list()
  for (l in 1:M){
    FF1[[l]]<- (t(y - b_l1[[l]]) %*% X_lsampled[[l]] - (t(t2) %*% t(X_lsampled[[l]]) %*% X_lsampled[[l]]))^2
  }
  FF2<- (1/(M * t1[3])) * (Reduce('+', FF1))
  
  ## Monte Carlo Approximation for the standard error of betag
  ## Since the betagsqu is block diagonal
  varBeta<- betagsqu[2,2] + FF2[,2] ## Observed Fisher Information, 
  standard_error <- sqrt(1/varBeta)  ##variance is given as inverse of the observed fisher information
  estimated_parameter <- t1[2]
  
  ### Wald test statistic
  compute_wald_statistic <- function(est_parameter, standard_error) {
    wald_statistic <- est_parameter / standard_error
    return(wald_statistic)
  }
  wald_value <- compute_wald_statistic(estimated_parameter, standard_error)
  degrees_of_freedom <- 1
  
  # Calculate p-value# pchisq by default compute left tail probabilities that's why we do
  # 1-p to obtain the corresponding right tail
  #p_value <- 1 - pchisq(wald_value, df = degrees_of_freedom)
  p_value <- pchisq(wald_value, df = degrees_of_freedom, lower.tail = FALSE)
  
  # Compare p-value to significance level (e.g., alpha = 0.05)
  #alpha <- 0.05
  # if (p_value < alpha) {
  #  print("Reject the null hypothesis: Parameter is significant")
  # } else {
  #    print("Fail to reject the null hypothesis: Parameter is not significant")
  #  }
  print(paste("p-value:", p_value))
  
  thetap <- c(t1,  p_value)
  
  #################################################################################
  #theta.hist<- as.data.frame(thetap)
  #print(theta.hist)
  #colnames(theta.hist) <- c("beta0", "beta1", "sigsqub", "sigsque", "p-value")
  #rownames(theta.hist) <- paste("it", 0:50, sep = "")
  return(list(theta.hat=thetahat, theta.hist = thetap))
}


#---#---#---# Call the EM function and write estimates to theta.hat.txt #---#---#---#---
#---#---#---# Read in data and Kinship matrix for testing:    #---#---#---#---#---#
#---#---#---# Make sure files are in the same folder to run the program in the cluster

##Test with single files:
#setwd("/Users/maryamonifade/Documents/EM_algorithm_files")
#Test_data2<- read.table("dat_comp1.txt", header = TRUE)
#kin_mat<- matrix(scan("kin1.txt"), nrow=5000, byrow=TRUE)
#theta1<- read.table("MYthetahat1.txt")
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

total_files=seq(from=1, to=1001, by=10) 

starting = total_files[task_id]
#Compute starting index

ending = total_files[task_id +1] - 1
#Compute ending index

for (it in (starting:ending)){
  
  #for (it in (101:177)){
  
  ### Read in the data
  dat_comp<- paste("/global/scratch/hpc4298/EM_EPS_data/GASTON_EM_files/dat_comp", it, ".txt", sep="")
  Test_data2<- read.table(dat_comp, header = TRUE)
  ###Use kinship matrices from gaston which has been internally standardised
  cov_matrice <- paste("/global/scratch/hpc4298/EM_EPS_data/GASTON_EM_files/kin", it, ".txt", sep="")
  kin_mat<-matrix(scan(cov_matrice), nrow=5000, byrow=TRUE)
  Em_estimate <- paste("MYthetahat", it,  ".txt", sep="")
  theta1<- read.table(Em_estimate, header=FALSE)
  
  ##EM-Wald  Model:
  Wald_test <- Em_Wald_test(Test_data2, M=100, kin_mat, theta1=theta1)
  modelEM.out<- Wald_test$theta.hist
  #write.table(modelEM.out, file =  paste("EMParams", it, ".txt", sep="") , col= T, row= FALSE, sep= "\t", quote=FALSE)
  write(modelEM.out, file = paste("P_value", it, ".txt", sep=""), sep= " ", append=FALSE, ncolumns=length(modelEM.out))
}





