### R script for Wald test for the randomly missing datasets
### January 11th

library(dplyr)

Em_Rand<- function(data, M=100, K, theta1){
  
  #First obtain genotype frequencies: 
  geno_comp<- data[,3]
  Xstar<- as.matrix(cbind(x0=1, x1=geno_comp)) # Complete genotype vector with intercept
  genotype_counts<- as.data.frame(table(geno_comp))  
  geno_prop<- c()
  for (i in genotype_counts["Freq"]){
    d<-length(geno_comp)
    geno_prop<- i/d
  }
  ###Initial estimates from GASTON
  theta0<- c()
  theta0[1:ncol(Xstar)]<- c(-0.0523, 0.0477) 
  theta0[3]<- 0.2  ###sigmasqub
  theta0[4]<- 1.001904 ##sigmasque
  theta0[5:7] <- geno_prop
  id <- data[,1]
  u.id <- unique(id)
  n <- length(u.id)
  y<- data[,2]
  beta0<- theta0[1:ncol(Xstar)]
  sigmasqub0<- theta0[3] ##Variance for the random effect
  sigmasque0<- theta0[4] ## Variance estimates of the residual error
  initial.theta <- theta0
  theta.hist <- theta0
  
  #for (it in 1:niter) {
    D <- sigmasqub0 * K
    invD<- Matrix::chol2inv(D)
    ## Variance of y
    idenRand<- diag(sigmasque0, nrow=n, ncol=n)
    Vary<- D + idenRand
    ## variance due to the suffcient statistics: check again
    Vy<- Matrix::chol2inv(invD + diag(rep(sigmasque0, n)))
    ## Variance of y
    #Vary<- D + diag(sigmasque0, nrow=n, ncol=n) 
    ## variance due to the sufficient statistics: check again
    #Vy<- solve(solve(D) + diag(rep(sigmasque0, n))) 
    #Vy<- solve(solve(D) + sigmasque0) 
    
    #-----||-----||-----||-----||-----||-----||-----||-----||-----
    #STEP2: Compute weights for the individuals with missing genotype data:
    
    ## For randomly missing individuals, we first randomly delete these observations
    ## maybe 5%, 10% or 20 % from the complete dataset.
    ## Then we obtain the IDS of these individuals so that we can insert the weights
    ## we obtain back into the full dataset after computing. 
    
    #Step 2(i): Randomly delete 5% of the genotypes: this is equal to about 125 observations
    ## The results in a smaller dataset now with 2407 observations.  
    ## This is saved in rand_data.
    #rand_data<- data[-sample(1:nrow(data), 125), ]
    
    #Adjust this line according to percentage of desired missingness
    rand_data<- data[-sample(1:nrow(data), (0.05 * nrow(data))), ]
    
    # Randomly sample observations to be deleted. 
    #samp_id<- sample(1:nrow(data), 125) ##IDs of obs to be deleted.
    
    # Randomly sample observations to be deleted. 
    samp_id<- sample(1:nrow(data), (0.05 * nrow(data))) ##IDs of obs to be deleted.
    
    ## Now use these ids to only delete from the genotype column
    obs_geno<- data$Geno[-(samp_id)] ##observed genotype
    ## The full observations missing are:
    dat_missing <- anti_join(data, rand_data, by=c("ID"))
    #Test that deleted genotypes are equal to the geno column in dat_missing. 
    ###setequal(obs_geno,  dat_missing$Geno) = TRUE
    
    ##Now we compute weights for individuals with missing genotypes
    ## Using the phenotype information for the individuals with missing genotypes:
    
    missing_geno_weights<- matrix(rep(0, nrow(dat_missing) * 3), ncol = 3)
    #Ind_var<- Vary[row(Vy) == col(Vy)]
    Ind_var<- Vary[row(Vary) == col(Vary)] ##Obtain the diagonal elements of the variance matrix 
    # Subset to obtain those for just the missing individuals
    Missing_ind_var<- Ind_var[as.integer(rownames(dat_missing))] 
    y_mis<- dat_missing$BMI
    
    miss_geno_probs<- matrix(rep(0, nrow(dat_missing) * 3), ncol = 3)
    for (i in seq_along(y_mis)) {
      for (j in 0:2){
        #y_imis <- y_mis[id_mis == u.idMis[i]]
        miss_geno_probs[i, j+1]<- dnorm(y_mis[i], (theta0[1] + j*theta0[2]), sqrt(Missing_ind_var[i])) * geno_prop[j+1]
      }
    }
    
    missing_geno_weights<- matrix(rep(0, nrow(dat_missing) * 3), ncol = 3)
    for (i in seq_along(y_mis)){
      for (j in seq_along(1:3)){
        missing_geno_weights[i,j]<- miss_geno_probs[i,j]/sum(miss_geno_probs[i,])
      }
    }
    
    #-----||-----||-----||-----||-----||-----||-----||-----||-----
    #STEP2b: Prepare the weights for the individuals with complete genotype data 
    # Recall the weight is 1 for the actual genotype of the individual with complete 
    # genotypes and 0 for the other two genotypes. We show this below:
    #-----||-----||-----||-----||-----||-----||-----||-----||----
    #First initialise 
    Weights_complete<- matrix(rep(0, nrow(rand_data) * 3), ncol = 3)
    for (i in 1:nrow(rand_data)){
      if (rand_data$Geno[i] =='0'){
        Weights_complete[i,1]<- 1
        Weights_complete[i,2]<- 0
        Weights_complete[i,3]<- 0
      } else if (rand_data$Geno[i] =='1'){
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
      for (l in seq(1,M)){
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
    
    g_obs<- c(rand_data$Geno)
    gl <-list()
    for (i in 1:nrow(sample_outcomes)){
      gl[[i]] <- c(g_obs, sample_outcomes[i,])
    }
    
    #-----//----------//----------//--------
    ## Now for each lth sampled genotype vector, combine it with the intercept column
    X_lsampled<- list()
    for (i in 1:nrow(sample_outcomes)){
      X_lsampled[[i]] <- as.matrix(cbind(Xstar[,1], gl[[i]]))
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
    XX<-solve(Reduce('+', D)) ###Taking the inverse of the sum of lists
    
    F1<- list()
    F2<- list()
    for (l in 1:M){
      F1[[l]]<- t(y - b_l1[[l]])
      F2[[l]] <- F1[[l]] %*% X_lsampled[[l]]
    }
    XZ<- t(Reduce('+', F2))
    betastar1<- ((XX) %*% XZ)
    #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
    # STEP2g:Estimating covariance parameters at the t+1th iteration.
    # At the t=0 iteration: 
    # For sigmaque
    #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
    #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
    
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
    
    ######First computation of this parameter (gives 34580).
    ###checked through again: not correct
    # #E<-0
    # ee <- c()
    # e<- list()
    # #ee<- list()
    # for (l in 1:M){
    #   e[[l]] <- (y - (X_lsampled[[l]] %*% betastar1) - b_l1[[l]])
    #   #ee[[l]]<- t(e[[l]]) %*% e[[l]]
    #   ee<- sum(t(e[[l]]) %*% e[[l]]) ##this is not really the sum of the list, only giving the last product
    # }
    # 
    # # E<- sum(ee)
    # #E<- Reduce('+', ee)
    # 
    # sigmasque1<-(sum(diag(Vy))) + ((1/M) * ee)
    
    #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
    # STEP2f:Estimating covariance parameters at the t+1th iteration.
    # At the t=0 iteration: 
    # For sigmasqub
    ##-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
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
    
    All_ind_weights<- as.data.frame(rbind(Weights_complete, missing_geno_weights))
    ## Now to compute the estimates of the genotype parameters:
    Obs_freq1<- length(which(All_ind_weights$V1 =="1"))
    miss_freq1<- All_ind_weights$V1[All_ind_weights$V1 < 1]
    gamma0<-  (Obs_freq1/n) + (sum(miss_freq1)/n)
    gamma1<- (length(which(All_ind_weights$V2 =="1"))/n) + sum(All_ind_weights$V2[All_ind_weights$V2 <1])/n
    gamma2<- 1-gamma0-gamma1
    geno_prop1 <- c(gamma0,gamma1,gamma2)
    
    ######## Collate all estimates #####
    thetahat <- c(betastar1,sigmasqub1,sigmasque1, geno_prop1)
    #theta.hist <- rbind(theta.hist, theta1)
    
    ### Monte carlo Appraoch to the Wald test.
    ### Second derivative for betag
    ## Gives a two by two matrix: To resuse this in the information matrix,
    ## we pick the diagonal elements which correspond to the second derivatives of 
    ## the intercept and betag
    t1<- c(theta1[,1], theta1[,2], theta1[,3], theta1[,4])
    t2 <- t1[1:2]
    
    Xll<- Reduce('+', D) ##
    betagsqu<- (1/(M * t1[3])) * Xll 
    
    ### First derivative computation:
    ##Initialize list to store components of the bracket
    FF1<- list()
    for (l in 1:M){
      FF1[[l]]<- (t(y - b_l1[[l]]) %*% X_lsampled[[l]] - (t(t2) %*% t(X_lsampled[[l]]) %*% X_lsampled[[l]]))^2
    }
    FF2<- (1/(M * t1[3])) * (Reduce('+', FF1))
    
    ## Monte Carlo Approximation for the standard error of betag
    ## Since the betagsqu is block diagonal 
    varBeta<- betagsqu[2,2] + FF2[,2]
    standard_error <- sqrt(1/varBeta)
    estimated_parameter <- t1[2]
    
    ### Wald test statistic
    compute_wald_statistic <- function(est_parameter, standard_error) {
     wald_statistic <- est_parameter / standard_error
     return(wald_statistic)
    }
    wald_value <- compute_wald_statistic(estimated_parameter, standard_error)
    #print(paste("Wald Test Statistic:", wald_value))
    
    # Obtain pvalue for the wald test statistic
    # Computed Wald test statistic
    #wald_value <- 4.27  # Replace this with your computed Wald test statistic
    
    # Degrees of freedom (df)
    degrees_of_freedom <- 1
    
    # Calculate p-value# pchisq by default compute left tail probabilities that's why we do 
    # 1-p to obtain the corresponding right tail
    #p_value <- 1 - pchisq(wald_value, df = degrees_of_freedom)
    p_value1 <- pchisq(wald_value, df = degrees_of_freedom, lower.tail = FALSE)
    
    # Compare p-value to significance level (e.g., alpha = 0.05)
    #alpha <- 0.05
    #if (p_value1 < alpha) {
    # print("Reject the null hypothesis: Parameter is significant")
    #} else {
    # print("Fail to reject the null hypothesis: Parameter is not significant")
    #}
    print(paste("p-value:", p_value1))
    thetap <- c(t1,  p_value1)
    
    return(list(theta.hat=thetahat, theta.hist = thetap))
    
  }

#---#---#---# Call the EM_wald_Randfunction #---#---#---#---


  slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
  task_id <- as.numeric(slurm_arrayid)

  total_files=seq(from=1, to=1001, by=10)  ###Total files here is 200 hence 69 arrays.

  starting = total_files[task_id]
  #Compute starting index

 ending = total_files[task_id +1] - 1
 #Compute ending index

 for (it in (starting:ending)){

#for (it in (101:1000)){
  
  ### Read in the data
  dat_comp<- paste("/global/scratch/hpc4298/EM_EPS_data/GASTON_EM_files/dat_comp", it, ".txt", sep="")
  Test_data2<- read.table(dat_comp, header = TRUE)
  ###Use kinship matrices from gaston which has been internally standardised
  cov_matrice <- paste("/global/scratch/hpc4298/EM_EPS_data/GASTON_EM_files/kin", it, ".txt", sep="")
  kin_mat<-matrix(scan(cov_matrice), nrow=5000, byrow=TRUE)
  Em_estimate <- paste("theta_rand", it,  ".txt", sep="")
  theta1<- read.table(Em_estimate, header=FALSE)
  
  ##EM Model:
  Wald_test<- Em_Rand(Test_data2, M=100, kin_mat, theta1)
  modelEM.out<- Wald_test$theta.hist
  write(modelEM.out, file = paste("P_vals", it, ".txt", sep=""), sep= " ", append=FALSE, ncolumns=length(modelEM.out))

}

 
 
 