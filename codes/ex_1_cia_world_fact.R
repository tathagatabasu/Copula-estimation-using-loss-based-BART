# packages
library(dplyr)
library(readr)

# dataset
cia_wf_data <- read.csv("countries.csv")

cia_wf_data <- cia_wf_data %>% dplyr::select(all_of(c("Country","People.and.Society..Life.expectancy.at.birth...male",
                                                      "People.and.Society..Life.expectancy.at.birth...female",
                                                      "People.and.Society..Literacy...male",
                                                      "People.and.Society..Literacy...female", 
                                                      "Economy..Real.GDP.per.capita")))

colnames(cia_wf_data) <- c("Country",
                              "Life_expectancy_M",
                              "Life_expectancy_F",
                              "Liter_M",
                              "Liter_F",
                              "GDP_PPP")

cia_wf_data = cia_wf_data %>%
  mutate(across(-c(Country),.fns = parse_number))

cia_wf_data <- na.omit(cia_wf_data)


U1_LE = ecdf(cia_wf_data$Life_expectancy_F)(cia_wf_data$Life_expectancy_F)
U2_LE = ecdf(cia_wf_data$Life_expectancy_M)(cia_wf_data$Life_expectancy_M)

U1_LT = ecdf(cia_wf_data$Liter_F)(cia_wf_data$Liter_F)
U2_LT = ecdf(cia_wf_data$Liter_M)(cia_wf_data$Liter_M)

par(mar = c(5,5,2,1), mfrow = c(1,3))

plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Life_expectancy_M, xlab = "GDP (log-scale)", ylab = "Life Expectancy (M)")
plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Life_expectancy_F, xlab = "GDP (log-scale)", ylab = "Life Expectancy (F)")
plot(U1_LE,U2_LE, xlab = "Life Expectancy (F)", ylab = "Life Expectancy (M)")

# 10, 3.5
par(mar = c(5,5,2,1), mfrow = c(1,3))

plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Liter_M, xlab = "GDP (log-scale)", ylab = "Literacy (M)")
plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Liter_F, xlab = "GDP (log-scale)", ylab = "Literacy (F)")
plot(U1_LT,U2_LT, xlab = "Literacy (F)", ylab = "Literacy (M)")

GDP <- as.data.frame((log(cia_wf_data$GDP_PPP) - min(log(cia_wf_data$GDP_PPP)))/(max(log(cia_wf_data$GDP_PPP)) - min(log(cia_wf_data$GDP_PPP))))
GDP <- as.matrix(GDP)
rownames(GDP) <- 1:nrow(GDP)

n.chain_par <- 10
n.iter_par <- 15000
n.born.out.par <- 1000
n.thin <- 1
incl.split_par <- TRUE
cont.unif_par <- TRUE
moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)

n.tree <- 10

################################################################################
# source files
################################################################################

source('MCMC_BART_copula.R')
source('import_functions.R')

lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
##########################################################

if(F){
  gauss_GDP_LE <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                         n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
  
  save(gauss_GDP_LE, file = paste0("gauss_gdp_LE_tree_",n.tree,".Rdata"))
  rm(gauss_GDP_LE)
  gc()
  
  t_GDP_LE <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                     n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
  
  save(t_GDP_LE, file = paste0("t_gdp_LE_tree_",n.tree,".Rdata"))
  rm(t_GDP_LE)
  gc()
  
  clayton_GDP_LE <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                           n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
                                           X = GDP,
                                           U1 = U1_LE,
                                           U2 = U2_LE,
                                           prior_list = lb.prior.def, 
                                           moves.prob = moves.prob_par, 
                                           starting.tree = NULL,
                                           cont.unif = cont.unif_par,
                                           include.split = incl.split_par,
                                           prop_mu = 0, prop_sigma = 0.2,
                                           theta_param_1 = 0, theta_param_2 = 1,
                                           var_param_1 = 1, var_param_2 = 2,
                                           prior_type = "N",
                                           cop_type = "clayton",
                                           adapt = F)
  
  save(clayton_GDP_LE, file = paste0("clayton_gdp_LE_tree_",n.tree,".Rdata"))
  rm(clayton_GDP_LE)
  gc()
  
  gumbel_GDP_LE <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                          n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
                                          cop_type = "gumbel",
                                          adapt = F)
  
  save(gumbel_GDP_LE, file = paste0("gumbel_gdp_LE_tree_",n.tree,".Rdata"))
  rm(gumbel_GDP_LE)
  gc()
  
  frank_GDP_LE <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                         n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
                                         cop_type = "frank",
                                         adapt = F)
  
  save(frank_GDP_LE, file = paste0("frank_gdp_LE_tree_",n.tree,".Rdata"))
  rm(frank_GDP_LE)
  gc()
}

if(T){
  gauss_GDP_LT <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                         n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
  
  save(gauss_GDP_LT, file = paste0("gauss_gdp_LT_tree_",n.tree,".Rdata"))
  rm(gauss_GDP_LT)
  gc()
  
  t_GDP_LT <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                     n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
  
  save(t_GDP_LT, file = paste0("t_gdp_LT_tree_",n.tree,".Rdata"))
  rm(t_GDP_LT)
  gc()
  
  clayton_GDP_LT <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                           n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
                                           cop_type = "clayton",
                                           adapt = F)
  
  save(clayton_GDP_LT, file = paste0("clayton_gdp_LT_tree_",n.tree,".Rdata"))
  rm(clayton_GDP_LT)
  gc()
  
  gumbel_GDP_LT <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                          n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
                                          cop_type = "gumbel",
                                          adapt = F)
  
  save(gumbel_GDP_LT, file = paste0("gumbel_gdp_LT_tree_",n.tree,".Rdata"))
  rm(gumbel_GDP_LT)
  gc()
  
  frank_GDP_LT <- multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                         n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
                                         cop_type = "frank",
                                         adapt = F)
  
  save(frank_GDP_LT, file = paste0("frank_gdp_LT_tree_",n.tree,".Rdata"))
  rm(frank_GDP_LT)
  gc()
  
}
