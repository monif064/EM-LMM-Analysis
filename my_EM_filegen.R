### Program for simulating files for running GEMMA programs  #####
### In this version, we include the candidate SNP into the   #####
### whole genome and just focus on that when obtaining the result.###
##Update: Modify program to just 

setwd("/global/home/hpc4298/GEMMA-0.98/myexample/GEMMA_quan_3")

#GEMMA_sim<- function(N, nSNP){

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

total_files=seq(from=1, to= 1001, by=10)

starting = total_files[task_id]
#Compute starting index

ending = total_files[task_id +1] - 1
#Compute ending index

for (k in (starting:ending)){
  
  numpop=2
  N=5000    #Total number of individuals
  nSNP=5000  #number of  snps
  Fst=0.01
  omega=c(0.5,0.5) ##MAF
  propnExtreme=0.1
  #nsim=100
  Fst.obs=vector(length=nSNP)
  pdiffs=vector(length=nSNP)
  genomat=matrix(nrow=N,ncol=nSNP)
  ##Simulate snps for each population
  for (i in 1:nSNP){
    p=runif(1,0.1,0.9) ##generating allele frequency p from uniform distribution
    alpha=p*(1-Fst)/Fst
    beta=(1-p)*(1-Fst)/Fst
    ps=rbeta(numpop,shape1=alpha,shape2=beta)
    
    for (j in 1:numpop){
      ind1=(j-1)*N*omega[j]+1 ##index values for individuals in population 1
      ind2=j*N*omega[j]  #individuals in population 2
      freqs=c(ps[j]^2,2*ps[j]*(1-ps[j]),(1-ps[j])^2)
      genomat[ind1:ind2,i]=sample(c(0,1,2),size=N*omega[j],replace=TRUE,prob=freqs)
    }
    
  }
  X<- genomat
  ####Candidate SNP generation
  geno_SNP<- matrix(nrow=N,ncol=1)
  p<- c(0.25,0.85)
  for (j in 1:numpop){
    ind1=(j-1)*N*omega[j]+1 ##index values for individuals in population 1
    ind2=j*N*omega[j]  #individuals in population 2
    freqs=c(p[j]^2,2*p[j]*(1-p[j]),(1-p[j])^2)
    geno_SNP[ind1:ind2,1]=sample(c(0,1,2),size=N*omega[j],replace=TRUE,prob=freqs)
  }
  ###Calculate allele frequencies of the CandSNP
  genotype_counts<- as.data.frame(table(geno_SNP))  
  geno_prop<- c()
  for (i in genotype_counts["Freq"]){
    d<-length(geno_SNP)
    geno_prop<- i/d
  }
  
 ##simulate phenotype values from normal distribution.
  ##simulate phenotypes dependent on population since we are testing for type 1 error. 
  pheno_indep <-c()
  pheno1<- c(pheno_indep, rnorm((ind1-1), mean= 0.07, sd=1))
  pheno2<- c(pheno_indep, rnorm((ind1:ind2), mean= -0.07, sd=1))
  pheno_indep<- c(pheno1,pheno2)
  ##individual ID's and 
  IND<- 1:N
  combined_indep <- cbind(IND, pheno_indep, geno_SNP, X) ### Complete dataset including the SNP for association 
  sorted_combined <- combined_indep[order(combined_indep[,2]),] ##sort to subset for EPS data
  
  ##Complete file for phenotypes, single SNP(not missing) and IDs. 
  ## i.e before subsetting for EPS data:write this out as well
  dat_complete<- sorted_combined[, 1:3]
  write.table(dat_complete, file = paste("dat_comp", k, ".txt", sep=""), col.names = FALSE, row.names = FALSE, sep= "\t", quote=FALSE)
  
  ## Obtain the EPSdata for the missing genotypes
  K = propnExtreme #proportion in the extremes
  Nums = nrow(sorted_combined)
  keep <- c(1:(K*Nums), (Nums-(K*Nums)+1):Nums)
  obs_data<- as.data.frame(sorted_combined[keep,])

  ### EPS data that we need in the second phase; Do not include the big genotype matrix
  ## Write this out to a file
  EPS_data<-obs_data[,c(1:3)]
  write.table(EPS_data, file = paste("eps_dat", k, ".txt", sep=""), col.names = FALSE, row.names = FALSE, sep= "\t", quote=FALSE)
  
  #epsdat<- c(rep(1,K*Nums),rep(0,K*Nums))
  #EPS_pheno <- as.matrix(cbind(epsdat,sorted_combined[keep,])) 
  
  ##isolate the genotype matrix corresponding to the EPS sampples and add the candidate SNP as the first SNP
  ##Not anymore: the genotype matrix here doesnt include the candidate SNP
  #EPS_geno<- EPS_pheno[,-c(1:3)]  
  #EPS_geno_cand <- rbind(geno_SNP, EPS_geno)
  #my_eps<- t(EPS_geno)
  
  ## Genotype matrix for GRM will be the X matrix showing all SNPs for the full individuals. 
  ## Just the X matrix above
  
  
  ### Below we try to create files for use in the GEMMA program to obtain a GRM matrix
  ##Naming the columns.
  pre<- "IND"; suf<- seq(1:ncol(X))
  colnames(X)<- paste(pre,suf, sep="")
  
  ##Naming the SNPs with rs names 
  prefix<- "rs"; suffix<- seq(from=1, to=nrow(X) , by = 1)
  SNP_IDs<- paste(prefix,suffix,sep="")
  ##alleles effect and reference 
  Ref<- c(rep("A", nrow(X)))
  Eff<- sample(c("C","G","T"), nrow(X), replace=TRUE, prob= NULL)
  EP<- cbind(SNP_IDs, Ref,Eff, X)
  write.table(EP, file = paste("genomat", k, ".geno.txt", sep=""), col.names = FALSE, row.names = FALSE, sep= "\t", quote=FALSE)
  
  
  ###Round the phenotypes to 4 decimal places.
  Pheno<- round(sorted_combined[,2], 4)
  
   ### Formatting a GEMMA type phenotype file
  ### A single line for the phenotypes: if binary, code cases 1 and controls 0
  #write.table(Pheno, file = "simu.pheno.txt", col.names = FALSE, row.names = FALSE, sep= "\n", quote=FALSE)
  write.table(Pheno, file = paste("genomat", k, ".pheno.txt", sep=""), col.names = FALSE, row.names = FALSE, sep= "\n", quote=FALSE)
  
  ###Formatting a GEMMA SNP annotation file
  #3 columns, SNP ID, basepair position, chromosome number
  
  POS_val= seq(1000, 5000, length.out = nrow(EP)) ##explain what lines does
  POS= format(POS_val,decimal.mark = ' ')  ##what does line do??
  SNP_anno<- cbind(SNP_IDs, POS, CHR=rep(1, nrow(EP)))
  #write.table(SNP_anno, file="simu.anno.txt", col.names = FALSE, row.names = FALSE, sep="\t", quote=FALSE)
  write.table(SNP_anno, file = paste("genomat", k, ".anno.txt", sep=""), col.names = FALSE, row.names = FALSE, sep= "\t", quote=FALSE)
  
  
  ##single candidate snp
  
  cand_ID<- c("rs243556"); r<- c("A"); E<- c("C")
  cand_gene<- c(cand_ID, r, E, geno_SNP)
  #cg<- as.matrix(cand_gene, nrow=1, ncol=length(geno_SNP)+3)
  write.table(t(cand_gene), file = paste("geno_cand", k, ".geno.txt", sep=""), col.names = FALSE, row.names = FALSE, sep= ",", quote=FALSE)
  
  #write.table(t(cand_gene), file="cand.geno.txt", col.names = FALSE, row.names=FALSE, sep=",", quote=FALSE)
  
}
