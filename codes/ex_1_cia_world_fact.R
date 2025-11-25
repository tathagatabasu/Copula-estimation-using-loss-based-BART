# packages
library(dplyr)
library(readr)
library(VineCopula)

set.seed(123)

# dataset
cia_wf_data <- read.csv("countries.csv")

cia_wf_data <- cia_wf_data %>% dplyr::select(all_of(c("Country","People.and.Society..Life.expectancy.at.birth...male",
                                                      "People.and.Society..Life.expectancy.at.birth...female",
                                                      "People.and.Society..Literacy...male",
                                                      "People.and.Society..Literacy...female",
                                                      "Economy..Real.GDP..purchasing.power.parity.")))
# "Economy..Real.GDP.per.capita")))

colnames(cia_wf_data) <- c("Country",
                           "Life_expectancy_M",
                           "Life_expectancy_F",
                           "Liter_M",
                           "Liter_F",
                           "GDP_PPP")

index_bil <- (grepl("billion ", cia_wf_data$GDP_PPP, fixed = TRUE))
index_tril <- (grepl("trillion ", cia_wf_data$GDP_PPP, fixed = TRUE))
index_mil <- (grepl("million ", cia_wf_data$GDP_PPP, fixed = TRUE))

cia_wf_data = cia_wf_data %>%
  mutate(across(-c(Country),.fns = parse_number))

cia_wf_data$GDP_PPP[index_bil] = 1000 * cia_wf_data$GDP_PPP[index_bil]
cia_wf_data$GDP_PPP[index_tril] = 1000000 * cia_wf_data$GDP_PPP[index_tril]

cia_wf_data <- na.omit(cia_wf_data)

plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Life_expectancy_M)
plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Life_expectancy_F)

plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Liter_M)
plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Liter_F)

U1_LE = pobs(cia_wf_data$Life_expectancy_F)
U2_LE = pobs(cia_wf_data$Life_expectancy_M)

U1_LT = pobs(cia_wf_data$Liter_F)
U2_LT = pobs(cia_wf_data$Liter_M)

plot(U1_LE,U2_LE)
plot(cia_wf_data$Life_expectancy_F,cia_wf_data$Life_expectancy_M)

plot(U1_LT,U2_LT)
plot(cia_wf_data$Liter_F,cia_wf_data$Liter_M)

GDP <- as.data.frame((log(cia_wf_data$GDP_PPP) - min(log(cia_wf_data$GDP_PPP)))/(max(log(cia_wf_data$GDP_PPP)) - min(log(cia_wf_data$GDP_PPP))))
GDP <- as.matrix(GDP)
rownames(GDP) <- 1:nrow(GDP)

n.chain_par <- 4
n.iter_par <- 30000
n.born.out.par <- 1000
n.thin <- 1
incl.split_par <- TRUE
cont.unif_par <- TRUE
moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)

n.tree <- 5

################################################################################
# source files
################################################################################

source('MCMC_BART_copula.R')
source('import_functions.R')

lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
##########################################################

if(T){
  gauss_GDP_LE_adapt <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                         n.tree = n.tree, n.chain = n.chain_par, n.cores = 5,
                                         X = GDP,
                                         U1 = U1_LE,
                                         U2 = U2_LE,
                                         prior_list = lb.prior.def, 
                                         moves.prob = moves.prob_par, 
                                         starting.tree = NULL,
                                         cont.unif = cont.unif_par,
                                         include.split = incl.split_par,
                                         prop_mu = 0, prop_sigma = .2,
                                         theta_param_1 = 0, theta_param_2 = 1,
                                         var_param_1 = 1, var_param_2 = 2,
                                         prior_type = "N",
                                         cop_type = "gauss",
                                         adapt = T)
  
  save(gauss_GDP_LE_adapt, file = paste0("gauss_GDP_LE_adapt_tree_",n.tree,".Rdata"))
  rm(gauss_GDP_LE_adapt)
  gc()
  gc()
  
  t_GDP_LE_adapt <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                     n.tree = n.tree, n.chain = n.chain_par, n.cores = 5,
                                     X = GDP,
                                     U1 = U1_LE,
                                     U2 = U2_LE,
                                     prior_list = lb.prior.def, 
                                     moves.prob = moves.prob_par, 
                                     starting.tree = NULL,
                                     cont.unif = cont.unif_par,
                                     include.split = incl.split_par,
                                     prop_mu = 0, prop_sigma = .2,
                                     theta_param_1 = 0, theta_param_2 = 1,
                                     var_param_1 = 1, var_param_2 = 2,
                                     prior_type = "N",
                                     cop_type = "t",
                                     adapt = T)
  
  save(t_GDP_LE_adapt, file = paste0("t_GDP_LE_adapt_tree_",n.tree,".Rdata"))
  rm(t_GDP_LE_adapt)
  gc()
  gc()
}

if(T){
  gauss_GDP_LT_adapt <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                         n.tree = n.tree, n.chain = n.chain_par, n.cores = 5,
                                         X = GDP,
                                         U1 = U1_LT,
                                         U2 = U2_LT,
                                         prior_list = lb.prior.def, 
                                         moves.prob = moves.prob_par, 
                                         starting.tree = NULL,
                                         cont.unif = cont.unif_par,
                                         include.split = incl.split_par,
                                         prop_mu = 0, prop_sigma = .2,
                                         theta_param_1 = 0, theta_param_2 = 1,
                                         var_param_1 = 1, var_param_2 = 2,
                                         prior_type = "N",
                                         cop_type = "gauss",
                                         adapt = T)
  
  save(gauss_GDP_LT_adapt, file = paste0("gauss_GDP_LT_adapt_tree_",n.tree,".Rdata"))
  rm(gauss_GDP_LT_adapt)
  gc()
  gc()
  
  t_GDP_LT_adapt <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                     n.tree = n.tree, n.chain = n.chain_par, n.cores = 5,
                                     X = GDP,
                                     U1 = U1_LT,
                                     U2 = U2_LT,
                                     prior_list = lb.prior.def, 
                                     moves.prob = moves.prob_par, 
                                     starting.tree = NULL,
                                     cont.unif = cont.unif_par,
                                     include.split = incl.split_par,
                                     prop_mu = 0, prop_sigma = .2,
                                     theta_param_1 = 0, theta_param_2 = 1,
                                     var_param_1 = 1, var_param_2 = 2,
                                     prior_type = "N",
                                     cop_type = "t",
                                     adapt = T)
  
  save(t_GDP_LT_adapt, file = paste0("t_GDP_LT_adapt_tree_",n.tree,".Rdata"))
  rm(t_GDP_LT_adapt)
  gc()
  gc()
}

if(T){
  gauss_GDP_LE <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                         n.tree = n.tree, n.chain = n.chain_par, n.cores = 5,
                                         X = GDP,
                                         U1 = U1_LE,
                                         U2 = U2_LE,
                                         prior_list = lb.prior.def, 
                                         moves.prob = moves.prob_par, 
                                         starting.tree = NULL,
                                         cont.unif = cont.unif_par,
                                         include.split = incl.split_par,
                                         prop_mu = 0, prop_sigma = .2,
                                         theta_param_1 = 0, theta_param_2 = 1,
                                         var_param_1 = 1, var_param_2 = 2,
                                         prior_type = "N",
                                         cop_type = "gauss",
                                         adapt = F)
  
  save(gauss_GDP_LE, file = paste0("gauss_GDP_LE_tree_",n.tree,".Rdata"))
  rm(gauss_GDP_LE)
  gc()
  gc()
  
  t_GDP_LE <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                     n.tree = n.tree, n.chain = n.chain_par, n.cores = 5,
                                     X = GDP,
                                     U1 = U1_LE,
                                     U2 = U2_LE,
                                     prior_list = lb.prior.def, 
                                     moves.prob = moves.prob_par, 
                                     starting.tree = NULL,
                                     cont.unif = cont.unif_par,
                                     include.split = incl.split_par,
                                     prop_mu = 0, prop_sigma = .2,
                                     theta_param_1 = 0, theta_param_2 = 1,
                                     var_param_1 = 1, var_param_2 = 2,
                                     prior_type = "N",
                                     cop_type = "t",
                                     adapt = F)
  
  save(t_GDP_LE, file = paste0("t_GDP_LE_tree_",n.tree,".Rdata"))
  rm(t_GDP_LE)
  gc()
  gc()
}

if(T){
  gauss_GDP_LT <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                         n.tree = n.tree, n.chain = n.chain_par, n.cores = 5,
                                         X = GDP,
                                         U1 = U1_LT,
                                         U2 = U2_LT,
                                         prior_list = lb.prior.def, 
                                         moves.prob = moves.prob_par, 
                                         starting.tree = NULL,
                                         cont.unif = cont.unif_par,
                                         include.split = incl.split_par,
                                         prop_mu = 0, prop_sigma = .2,
                                         theta_param_1 = 0, theta_param_2 = 1,
                                         var_param_1 = 1, var_param_2 = 2,
                                         prior_type = "N",
                                         cop_type = "gauss",
                                         adapt = F)
  
  save(gauss_GDP_LT, file = paste0("gauss_GDP_LT_tree_",n.tree,".Rdata"))
  rm(gauss_GDP_LT)
  gc()
  gc()
  
  t_GDP_LT <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                     n.tree = n.tree, n.chain = n.chain_par, n.cores = 5,
                                     X = GDP,
                                     U1 = U1_LT,
                                     U2 = U2_LT,
                                     prior_list = lb.prior.def, 
                                     moves.prob = moves.prob_par, 
                                     starting.tree = NULL,
                                     cont.unif = cont.unif_par,
                                     include.split = incl.split_par,
                                     prop_mu = 0, prop_sigma = .2,
                                     theta_param_1 = 0, theta_param_2 = 1,
                                     var_param_1 = 1, var_param_2 = 2,
                                     prior_type = "N",
                                     cop_type = "t",
                                     adapt = F)
  
  save(t_GDP_LT, file = paste0("t_GDP_LT_tree_",n.tree,".Rdata"))
  rm(t_GDP_LT)
  gc()
  gc()
}
