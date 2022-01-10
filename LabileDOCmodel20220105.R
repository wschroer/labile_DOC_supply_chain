#load required packages
library(truncnorm)

#create variable names vector
names <- {c( 
  'SoS1',
  'SoS2',
  'SoS_tot',
  'BP_r',
  'fract_recycled',
  'BP.PP2',
  'PER_r',
  'BGE_r',
  'lys_P_r', 
  'micro_graz_r',
  'slop',
  'zoo_DOC_egs',
  'zoo_AE',
  'micro_GGE_f_r',
  'micro_GGE_d_r',
  'micro_GGE_c_r',
  'micro_GGE_r',
  'transfer_eff_r',
  'meso_GGE_r',
  'exc',
  'PPP',
  'PPP_aval',
  'BP_r',
  'lys_B_TPP',
  'micro_ing_r',
  'micro_ass_r',
  'micro_prod_r',
  'micro_egs_r',
  'micro_exc_r',
  'micro_resp_r',
  'meso_grazed_r',
  'meso_slop_r',
  'meso_ing_r',
  'meso_ass_r',
  'meso_prod_r',
  'meso_egs_r',
  'meso_exc_r',
  'meso_resp_r',
  'lys_B_r',
  'micro_ing2_r',
  'micro_ass2_r',
  'micro_prod2_r',
  'micro_egs2_r',
  'micro_exc2_r',
  'micro_resp2_r',
  'meso_grazed2_r',
  'meso_slop2_r',
  'meso_ing2_r',
  'meso_ass2_r',
  'meso_prod2_r',
  'meso_egs2_r',
  'meso_exc2_r',
  'meso_resp2_r',
  'micro_ing_tot',
  'micro_ass_tot',
  'micro_prod_tot',
  'micro_egs_tot',
  'micro_exc_tot',
  'micro_resp_tot',
  'meso_grazed_tot',
  'meso_slop_tot',
  'meso_ing_tot',
  'meso_ass_tot',
  'meso_prod_tot',
  'meso_egs_tot',
  'meso_exc_tot',
  'meso_resp_tot'
)}

#Define input parameters

TPP <- 1 #total primary production (TPP)
PER <- (read.csv('PER.csv')[ ,1])/100 #percent extracellular release as a fraction of TPP (PER), values compiled from literature review
BGE <- read.csv('BGE.csv')[ ,1] #bacterial growth efficiency (BGE), values compiled from literature review
lys_P <- c(0.02, 0.1) #fraction of TPP released by phage lysis of phytoplankton, range from Wilhelm and Suttle 1999
micro_graz <- read.csv('microGraz.csv', header=TRUE, sep=',')[ ,2]/100 #fraction of particulate primary production (PPP) grazed by microzooplankton, values from Schmoker et al. 2013
zoo_DOC_egs <- 0.5 #fraction of zooplankton egestia that is DOC from Urban-Rich 1999
zoo_AE <- 0.7 #zooplankton assimilation efficiency from Steinberg and Landry 2017
micro_GGE_flagel <- c(0.32, 0.17, 0.10, 0.63) #values describing the distribution (mean, sd, min, max) of microzooplankton (flagellate) gross growth efficiencies (GGE) calculated from data in Straile 1997
micro_GGE_dino <- c(0.30, 0.21, 0.04, 0.67) #same as above for dinoflagellates
micro_GGE_cili <- c(0.26, 0.12, 0.09, 0.48) #same as above for ciliates
meso_GGE <- c(0.26, 0.21, 0.01, 0.68) #same as above for mesozooplankton
max_slop <- 0.296 #maximum value of mesozooplankton sloppy feeding from Møller 2007
exc <- c(0.12) #fraction of carbon ingested by zooplankton that is excreted, value from Saba et al. 2011b
lys_B <- c(0.2, 0.3) #fraction of bacterial carbon lysed by viruses, range max and min from Fuhrman 1999; Wilhelm and Suttle 1999

set.seed(1421)
out <- data.frame(matrix(data=NaN, ncol=length(names))) #define data frame for model output 
colnames(out) <- names

#model code
# '_r' suffix on parameters indicates they are derived from the random draw
n_runs <- 1000 #number of model runs, for code testing set to 1,000 
x <- c()
for(i in 1:n_runs){ 
  #Primary parameters
  PER_r <- sample(PER, 1, replace=TRUE)
  BGE_r <- sample(BGE, 1, replace=TRUE)
  lys_P_r <- runif(1,lys_P[1],lys_P[2])
  micro_graz_r <- sample(micro_graz, 1, replace=TRUE)
  slop <- runif(1, 0, max_slop)
  micro_GGE_f_r <- rtruncnorm(1, a=micro_GGE_flagel[3], b=micro_GGE_flagel[4], mean=micro_GGE_flagel[1], sd=micro_GGE_flagel[2]) #draw from truncated normal distribution 
  micro_GGE_d_r <- rtruncnorm(1, a=micro_GGE_dino[3], b=micro_GGE_dino[4], mean=micro_GGE_dino[1], sd=micro_GGE_dino[2]) #draw from truncated normal distribution 
  micro_GGE_c_r <- rtruncnorm(1, a=micro_GGE_cili[3], b=micro_GGE_cili[4], mean=micro_GGE_cili[1], sd=micro_GGE_cili[2]) #draw from truncated normal distribution 
  micro_GGE_r <- mean(c(micro_GGE_f_r, micro_GGE_d_r, micro_GGE_c_r)) #mean microzooplankton GGE to use in subsequent calculations  
  transfer_eff_r <- runif(1, 0.3, 1) #microzooplanton C transfer efficincy to mesozooplankton from Steinberg and Landry 2017
  meso_GGE_r <- rtruncnorm(1, a=meso_GGE[3], b=meso_GGE[4], mean=meso_GGE[1], sd=meso_GGE[2]) #draw from truncated normal distribution 
  lys_B_r <- runif(1,lys_B[1],lys_B[2]) 
  
  #Derived Parameters, primary production circuit
  PPP <- TPP - PER_r
  PPP_aval <- TPP * PPP - lys_P_r
  micro_ing_r <- micro_graz_r * PPP_aval
  micro_ass_r <- micro_ing_r * zoo_AE
  micro_prod_r <- micro_ing_r * micro_GGE_r
  micro_egs_r <- (micro_ing_r - micro_ass_r) * zoo_DOC_egs
  micro_exc_r <- micro_ing_r * exc
  micro_resp_r <- micro_ass_r - micro_prod_r - micro_exc_r
  meso_grazed_r <- PPP_aval - (micro_graz_r * PPP_aval) + transfer_eff_r * micro_prod_r
  meso_slop_r <- meso_grazed_r * slop
  meso_ing_r <- meso_grazed_r - meso_slop_r
  meso_ass_r <- meso_ing_r * zoo_AE
  meso_prod_r <- meso_ing_r * meso_GGE_r
  meso_egs_r <- (meso_ing_r - meso_ass_r) * zoo_DOC_egs
  meso_exc_r <- meso_ing_r * exc 
  meso_resp_r <- meso_ass_r - meso_prod_r - meso_exc_r
  
  #calculate DOC sum of sources, primary circuit 
  SoS1 <- PER_r + lys_P_r + micro_egs_r + micro_exc_r + meso_slop_r + meso_egs_r + meso_exc_r
  
  #initate the recycled C circuit with bacterial production (BP) caclulated from primary sum of sources
  BP_r <- SoS1 * BGE_r 
  BP.PP2 <- BP_r / PPP
  lys_B_TPP <- BP_r * lys_B_r
  micro_ing2_r <- BP_r - lys_B_TPP
  micro_ass2_r <- micro_ing2_r * zoo_AE
  micro_prod2_r <- micro_ing2_r * micro_GGE_r
  micro_egs2_r <- (micro_ing2_r - micro_ass2_r) * zoo_DOC_egs
  micro_exc2_r <- micro_ing2_r * exc
  micro_resp2_r <- micro_ass2_r - micro_prod2_r - micro_exc2_r
  meso_grazed2_r <- transfer_eff_r * micro_prod2_r
  meso_slop2_r <- meso_grazed2_r * slop
  meso_ing2_r <- meso_grazed2_r - meso_slop2_r
  meso_ass2_r <- meso_ing2_r * zoo_AE
  meso_prod2_r <- meso_ing2_r * meso_GGE_r
  meso_egs2_r <- (meso_ing2_r - meso_ass2_r) * zoo_DOC_egs
  meso_exc2_r <- meso_ing2_r * exc
  meso_resp2_r <- meso_ass2_r - meso_prod2_r - meso_exc2_r
  
  #calculate DOC sum of sources from the recycled C circuit
  SoS2 <- lys_B_TPP +  micro_egs2_r + micro_exc2_r + meso_slop2_r + meso_egs2_r + meso_exc2_r
  
  #calculate total values as the sum of primary and recycled circuits 
  micro_ing_tot <- micro_ing_r + micro_ing2_r
  micro_ass_tot <- micro_ass_r + micro_ass2_r
  micro_prod_tot <- micro_prod_r + micro_prod2_r
  micro_egs_tot <- micro_egs_r + micro_egs2_r
  micro_exc_tot <- micro_exc_r + micro_exc2_r
  micro_resp_tot <- micro_resp_r + micro_resp2_r
  meso_grazed_tot <- meso_grazed_r + meso_grazed2_r
  meso_slop_tot <- meso_slop_r + meso_slop2_r
  meso_ing_tot <- meso_ing_r + meso_ing2_r
  meso_ass_tot <- meso_ass_r + meso_ass2_r
  meso_prod_tot <- meso_prod_r + meso_prod2_r
  meso_egs_tot <- meso_egs_r + meso_egs2_r
  meso_exc_tot <- meso_exc_r + meso_exc2_r
  meso_resp_tot <- meso_resp_r + meso_resp2_r
  SoS_tot <- SoS1 + SoS2
  fract_recycled <- SoS2 / SoS_tot
  
  for(j in 1:length(names)){
    x[j] <- eval(parse(text=names[j]))
  }
  out[i, ] <- x #save parameter and variable values for each run to output data frame
}

out <- out[out$meso_resp_r>0 & out$meso_resp2_r>0 & out$micro_resp_r>0 & out$micro_resp2_r>0, ] #filter out negative respiration values
write.csv(out, file='out_NAME_DATE.csv')

colSummary <- function(data){
  tab <- data[1:5, ]
  rownames(tab) <- c('mean', 'sd','median', 'min', 'max')
  tab[1, ] <- colMeans(data)
  for(i in 1:ncol(data)){
    tab[2, i] <- sd(data[ ,i])
  }
  for(i in 1:ncol(data)){
    tab[3, i] <- median(data[ ,i])
  }
  for(i in 1:ncol(data)){
    tab[4, i] <- min(data[ ,i])
  }
  for(i in 1:ncol(data)){
    tab[5, i] <- max(data[ ,i])
  }
  tab
}

out_summ <- colSummary(out)
out_summ_trans <- as.data.frame(t(as.matrix(out_summ)))
write.csv(out_summ_trans, file='model_summ_NAME_DATE.csv')
