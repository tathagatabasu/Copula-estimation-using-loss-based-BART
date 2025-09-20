# packages
library(data.tree)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(VineCopula)
library(MASS)   # For multivariate normal functions
library(coda)   # For MCMC diagnostics
library(plot3D)
library(gplots)
library(xtable)
require(foreach)
require(parallel)
require(doParallel)
library(calculus)

# dataset
cia_wf_data <- read.csv("countries.csv")

colnames_for_anlysis <- c("Country","People.and.Society..Life.expectancy.at.birth...male",
                          "People.and.Society..Life.expectancy.at.birth...female",
                          "People.and.Society..Literacy...male",
                          "People.and.Society..Literacy...female", 
                          "Economy..Real.GDP.per.capita")

cia_wf_data_le_vs_gdp <- cia_wf_data %>% dplyr::select(all_of(colnames_for_anlysis))

index_2023 <- (grepl("(2023 est.)", cia_wf_data_le_vs_gdp$Economy..Real.GDP.per.capita, fixed = TRUE))


cia_wf_data_2023 <- cia_wf_data_le_vs_gdp[index_2023,]

# index_bil <- (grepl("billion (2023 est.)", cia_wf_data_2023$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))
# index_tril <- (grepl("trillion (2023 est.)", cia_wf_data_2023$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))
# index_mil <- (grepl("million (2023 est.)", cia_wf_data_2023$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))

index <- (grepl(" (2023 est.)", cia_wf_data_2023$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))

colnames(cia_wf_data_2023) <- c("Country",
                              "Life_expectancy_M",
                              "Life_expectancy_F",
                              "Liter_M",
                              "Liter_F",
                              "GDP_PPP")

library(readr)

cia_wf_data_2023 = cia_wf_data_2023 %>%
  mutate(across(-c(Country),.fns = parse_number))

# commented as repating this will be problematic

# cia_wf_data_2023$GDP_PPP[index_bil] = 1000 * cia_wf_data_2023$GDP_PPP[index_bil]
# cia_wf_data_2023$GDP_PPP[index_tril] = 1000000 * cia_wf_data_2023$GDP_PPP[index_tril]

cia_wf_data_2023 <- na.omit(cia_wf_data_2023)

plot(log(cia_wf_data_2023$GDP_PPP),cia_wf_data_2023$Life_expectancy_M)
plot(log(cia_wf_data_2023$GDP_PPP),cia_wf_data_2023$Life_expectancy_F)

plot(log(cia_wf_data_2023$GDP_PPP),cia_wf_data_2023$Liter_M)
plot(log(cia_wf_data_2023$GDP_PPP),cia_wf_data_2023$Liter_F)

U1_LE = ecdf(cia_wf_data_2023$Life_expectancy_F)(cia_wf_data_2023$Life_expectancy_F)
U2_LE = ecdf(cia_wf_data_2023$Life_expectancy_M)(cia_wf_data_2023$Life_expectancy_M)

U1_LT = ecdf(cia_wf_data_2023$Liter_F)(cia_wf_data_2023$Liter_F)
U2_LT = ecdf(cia_wf_data_2023$Liter_M)(cia_wf_data_2023$Liter_M)

plot(U1_LE,U2_LE)
plot(cia_wf_data_2023$Life_expectancy_F,cia_wf_data_2023$Life_expectancy_M)

plot(U1_LT,U2_LT)
plot(cia_wf_data_2023$Liter_F,cia_wf_data_2023$Liter_M)

GDP <- as.data.frame((log(cia_wf_data_2023$GDP_PPP) - min(log(cia_wf_data_2023$GDP_PPP)))/(max(log(cia_wf_data_2023$GDP_PPP)) - min(log(cia_wf_data_2023$GDP_PPP))))
GDP <- as.matrix(GDP)
rownames(GDP) <- 1:nrow(GDP)

n.chain_par <- 10
n.iter_par <- 6000
n.born.out.par <- 1000
n.thin <- 1
incl.split_par <- TRUE
cont.unif_par <- TRUE
moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)

n.tree <- 10

################################################################################
# source files
################################################################################

# source('mclapply.R') # if run on windows uncomment it

# source('MCMC_BART_copula.R')
source('import_functions.R')
source('test_MCMC_copula_mult_tree.R')

lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
##########################################################

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
                                 prop_mu = 0, prop_sigma = .2/n.tree,
                                 theta_param_1 = 0, theta_param_2 = .3,
                                 var_param_1 = 1, var_param_2 = 2,
                                 prior_type = "N",
                                 cop_type = "gauss",
                                adapt = T)

save(gauss_GDP_LE, file = "gaussgdp_LE.Rdata")
rm(gauss_GDP_LE)

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
                            prop_mu = 0, prop_sigma = rep(.2/n.tree,n.tree),
                            theta_param_1 = 0, theta_param_2 = .3,
                            var_param_1 = 1, var_param_2 = 2,
                            prior_type = "N",
                            cop_type = "frank",
                            adapt = T)

save(frank_GDP_LE, file = "frankgdp_LE.Rdata")
rm(frank_GDP_LE)

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
                            prop_mu = 0, prop_sigma = rep(.2/n.tree,n.tree),
                            theta_param_1 = 0, theta_param_2 = .3,
                            var_param_1 = 1, var_param_2 = 2,
                            prior_type = "N",
                            cop_type = "t",
                            adapt = T)

save(t_GDP_LE, file = "tgdp_LE.Rdata")
rm(t_GDP_LE)

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
                                  prop_mu = 0, prop_sigma = rep(.2/n.tree,n.tree),
                                  theta_param_1 = 0, theta_param_2 = .3,
                                  var_param_1 = 1, var_param_2 = 2,
                                  prior_type = "N",
                                  cop_type = "clayton",
                                  adapt = T)

save(clayton_GDP_LE, file = "claytongdp_LE.Rdata")
rm(clayton_GDP_LE)

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
                                             prop_mu = 0, prop_sigma = rep(.2/n.tree,n.tree),
                                             theta_param_1 = 0, theta_param_2 = .3,
                                             var_param_1 = 1, var_param_2 = 2,
                                             prior_type = "N",
                                             cop_type = "gumbel",
                                             adapt = T)

save(gumbel_GDP_LE, file = "gumbelgdp_LE.Rdata")
rm(gumbel_GDP_LE)

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
                                       prop_mu = 0, prop_sigma = .2/n.tree,
                                       theta_param_1 = 0, theta_param_2 = .3,
                                       var_param_1 = 1, var_param_2 = 2,
                                       prior_type = "N",
                                       cop_type = "gauss",
                                       adapt = T)

save(gauss_GDP_LT, file = "gaussgdp_LT.Rdata")
rm(gauss_GDP_LT)

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
                                       prop_mu = 0, prop_sigma = rep(.2/n.tree,n.tree),
                                       theta_param_1 = 0, theta_param_2 = .3,
                                       var_param_1 = 1, var_param_2 = 2,
                                       prior_type = "N",
                                       cop_type = "frank",
                                       adapt = T)

save(frank_GDP_LT, file = "frankgdp_LT.Rdata")
rm(frank_GDP_LT)

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
                                   prop_mu = 0, prop_sigma = rep(.2/n.tree,n.tree),
                                   theta_param_1 = 0, theta_param_2 = .3,
                                   var_param_1 = 1, var_param_2 = 2,
                                   prior_type = "N",
                                   cop_type = "t",
                                   adapt = T)

save(t_GDP_LT, file = "tgdp_LT.Rdata")
rm(t_GDP_LT)

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
                                         prop_mu = 0, prop_sigma = rep(.2/n.tree,n.tree),
                                         theta_param_1 = 0, theta_param_2 = .3,
                                         var_param_1 = 1, var_param_2 = 2,
                                         prior_type = "N",
                                         cop_type = "clayton",
                                         adapt = T)

save(clayton_GDP_LT, file = "claytongdp_LT.Rdata")
rm(clayton_GDP_LT)

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
                                        prop_mu = 0, prop_sigma = rep(.2/n.tree,n.tree),
                                        theta_param_1 = 0, theta_param_2 = .3,
                                        var_param_1 = 1, var_param_2 = 2,
                                        prior_type = "N",
                                        cop_type = "gumbel",
                                        adapt = T)

save(gumbel_GDP_LT, file = "gumbelgdp_LT.Rdata")
rm(gumbel_GDP_LT)

################################################################################

if(F){
  gauss_list_pred_lb <- lapply(1:length(gauss_GDP_LE$trees), \(idx) BART_calculate_pred(gauss_GDP_LE$trees[[idx]], GDP))
  
  gauss_pred_val = do.call(rbind,gauss_list_pred_lb)
  
  gauss_pred_val_vec = as.vector(gauss_pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  frank_list_pred_lb <- lapply(1:length(frank_GDP_LE$trees), \(idx) BART_calculate_pred(frank_GDP_LE$trees[[idx]], GDP))
  
  frank_pred_val = do.call(rbind,frank_list_pred_lb)
  
  frank_pred_val_vec = as.vector(frank_pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  t_list_pred_lb <- lapply(1:length(t_GDP_LE$trees), \(idx) BART_calculate_pred(t_GDP_LE$trees[[idx]], GDP))
  
  t_pred_val = do.call(rbind,t_list_pred_lb)
  
  t_pred_val_vec = as.vector(t_pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  clayton_list_pred_lb <- lapply(1:length(clayton_GDP_LE$trees), \(idx) BART_calculate_pred(clayton_GDP_LE$trees[[idx]], GDP))
  
  clayton_pred_val = do.call(rbind,clayton_list_pred_lb)
  
  clayton_pred_val_vec = as.vector(clayton_pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  gumbel_list_pred_lb <- lapply(1:length(gumbel_GDP_LE$trees), \(idx) BART_calculate_pred(gumbel_GDP_LE$trees[[idx]], GDP))
  
  gumbel_pred_val = do.call(rbind,gumbel_list_pred_lb)
  
  gumbel_pred_val_vec = as.vector(gumbel_pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  # obs
  pred_obs = rep(GDP, each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  # like
  
  gauss_like_val <- apply(gauss_pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1_LE, U2_LE))
  
  gauss_like_df <-data.frame("nn" = gauss_like_val)
  gauss_like_df$idx <- rep(1:n.iter_par, n.chain_par)
  gauss_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  gauss_pl_like <- gauss_like_df %>%
    filter(idx > 2000) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  gauss_pl_like
  
  frank_like_val <- apply(frank_pred_val, 1, function(x)loglik_frank(link_frank(x), U1_LE, U2_LE))
  
  frank_like_df <-data.frame("nn" = frank_like_val)
  frank_like_df$idx <- rep(1:n.iter_par, n.chain_par)
  frank_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  frank_pl_like <- frank_like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  frank_pl_like
  
  t_like_val <- apply(t_pred_val, 1, function(x)loglik_t(link_t(x), U1_LE, U2_LE))
  
  t_like_df <-data.frame("nn" = t_like_val)
  t_like_df$idx <- rep(1:n.iter_par, n.chain_par)
  t_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  t_pl_like <- t_like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  t_pl_like
  
  clayton_like_val <- apply(clayton_pred_val, 1, function(x)loglik_clayton(link_clayton(x), U1_LE, U2_LE))
  
  clayton_like_df <-data.frame("nn" = clayton_like_val)
  clayton_like_df$idx <- 1:(n.chain_par*n.iter_par)
  
  clayton_like_df <-data.frame("nn" = clayton_like_val)
  clayton_like_df$idx <- rep(1:n.iter_par, n.chain_par)
  clayton_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  clayton_pl_like <- clayton_like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  clayton_pl_like
  
  gumbel_like_val <- apply(gumbel_pred_val, 1, function(x)loglik_gumbel(link_gumbel(x), U1_LE, U2_LE))
  
  gumbel_like_df <-data.frame("nn" = gumbel_like_val)
  gumbel_like_df$idx <- 1:(n.chain_par*n.iter_par)
  
  gumbel_like_df <-data.frame("nn" = gumbel_like_val)
  gumbel_like_df$idx <- rep(1:n.iter_par, n.chain_par)
  gumbel_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  gumbel_pl_like <- gumbel_like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  gumbel_pl_like
  
  # dataframe for plotting
  pred_cond <- data.frame("obs" = pred_obs)
  pred_cond$obs = pred_obs
  pred_cond$gauss_y = link_gauss(gauss_pred_val_vec)
  pred_cond$gauss_tau = BiCopPar2Tau(1, pred_cond$gauss_y)
  pred_cond$frank_y = link_frank(frank_pred_val_vec)
  pred_cond$frank_tau = BiCopPar2Tau(5, pred_cond$frank_y)
  pred_cond$t_y = link_t(t_pred_val_vec)
  pred_cond$t_tau = BiCopPar2Tau(2, pred_cond$t_y)
  pred_cond$clayton_y = link_clayton(clayton_pred_val_vec)
  pred_cond$clayton_tau = BiCopPar2Tau(3, pred_cond$clayton_y)
  pred_cond$gumbel_y = link_gumbel(gumbel_pred_val_vec)
  pred_cond$gumbel_tau = BiCopPar2Tau(4, pred_cond$gumbel_y)
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(obs) %>%
    summarise(gauss_theta_mean = mean(gauss_y), gauss_theta_q975 = quantile(gauss_y, .975), gauss_theta_q025 = quantile(gauss_y, .025),
              frank_theta_mean = mean(frank_y), frank_theta_q975 = quantile(frank_y, .975), frank_theta_q025 = quantile(frank_y, .025),
              t_theta_mean = mean(t_y), t_theta_q975 = quantile(t_y, .975), t_theta_q025 = quantile(t_y, .025),
              clayton_theta_mean = mean(clayton_y), clayton_theta_q975 = quantile(clayton_y, .975), clayton_theta_q025 = quantile(clayton_y, .025),
              gumbel_theta_mean = mean(gumbel_y), gumbel_theta_q975 = quantile(gumbel_y, .975), gumbel_theta_q025 = quantile(gumbel_y, .025)) 
  
  pred_cond_mod_tau = pred_cond_thin %>%
    group_by(obs) %>%
    summarise(gauss_tau_mean = mean(gauss_tau), gauss_tau_q975 = quantile(gauss_tau, .975), gauss_tau_q025 = quantile(gauss_tau, .025),
              frank_tau_mean = mean(frank_tau), frank_tau_q975 = quantile(frank_tau, .975), frank_tau_q025 = quantile(frank_tau, .025),
              t_tau_mean = mean(t_tau), t_tau_q975 = quantile(t_tau, .975), t_tau_q025 = quantile(t_tau, .025),
              clayton_tau_mean = mean(clayton_tau), clayton_tau_q975 = quantile(clayton_tau, .975), clayton_tau_q025 = quantile(clayton_tau, .025),
              gumbel_tau_mean = mean(gumbel_tau), gumbel_tau_q975 = quantile(gumbel_tau, .975), gumbel_tau_q025 = quantile(gumbel_tau, .025)) 
  
  ggplot(pred_cond_mod_tau) +
    geom_line(aes(obs, gauss_tau_mean, col = "gauss")) +
    geom_line(aes(obs, gauss_tau_q975, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, gauss_tau_q025, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, frank_tau_mean, col = "frank")) +
    geom_line(aes(obs, frank_tau_q975, col = "frank"),linetype="dotted") +
    geom_line(aes(obs, frank_tau_q025, col = "frank"),linetype="dotted") +
    geom_line(aes(obs, t_tau_mean, col = "t")) +
    geom_line(aes(obs, t_tau_q975, col = "t"),linetype="dotted") +
    geom_line(aes(obs, t_tau_q025, col = "t"),linetype="dotted") +
    geom_line(aes(obs, clayton_tau_mean, col = "clayton")) +
    geom_line(aes(obs, clayton_tau_q975, col = "clayton"),linetype="dotted") +
    geom_line(aes(obs, clayton_tau_q025, col = "clayton"),linetype="dotted") +
    geom_line(aes(obs, gumbel_tau_mean, col = "gumbel")) +
    geom_line(aes(obs, gumbel_tau_q975, col = "gumbel"),linetype="dotted") +
    geom_line(aes(obs, gumbel_tau_q025, col = "gumbel"),linetype="dotted") +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  cor(U1_LE,U2_LE, method = "kendall")
  cor(U1_LT,U2_LT, method = "kendall")
  
  GDP_gauss_pred <- BiCopSim(N = nrow(GDP), family = 1, par = link_gauss(colMeans(gauss_pred_val)))
  
  gauss_pred_U1_LE = GDP_gauss_pred[,1]
  gauss_pred_U2_LE = GDP_gauss_pred[,2]
  
  plot(U1_LE,U2_LE)
  plot(gauss_pred_U1_LE,gauss_pred_U2_LE)
  
  GDP_frank_pred <- BiCopSim(N = nrow(GDP), family = 5, par = link_frank(colMeans(frank_pred_val)))
  
  frank_pred_U1_LE = GDP_frank_pred[,1]
  frank_pred_U2_LE = GDP_frank_pred[,2]
  
  plot(U1_LE,U2_LE)
  plot(frank_pred_U1_LE,frank_pred_U2_LE)
  
  GDP_t_pred <- BiCopSim(N = nrow(GDP), family = 2, par = link_t(colMeans(t_pred_val)),par2 = 3)
  
  t_pred_U1_LE = GDP_t_pred[,1]
  t_pred_U2_LE = GDP_t_pred[,2]
  
  plot(U1_LE,U2_LE)
  plot(t_pred_U1_LE,t_pred_U2_LE)
  
  GDP_clayton_pred <- BiCopSim(N = nrow(GDP), family = 3, par = link_clayton(colMeans(clayton_pred_val)))
  
  clayton_pred_U1_LE = GDP_clayton_pred[,1]
  clayton_pred_U2_LE = GDP_clayton_pred[,2]
  
  plot(U1_LE,U2_LE)
  plot(clayton_pred_U1_LE,clayton_pred_U2_LE)
  
  GDP_gumbel_pred <- BiCopSim(N = nrow(GDP), family = 4, par = link_gumbel(colMeans(gumbel_pred_val)))
  
  gumbel_pred_U1_LE = GDP_gumbel_pred[,1]
  gumbel_pred_U2_LE = GDP_gumbel_pred[,2]
  
  plot(U1_LE,U2_LE)
  plot(gumbel_pred_U1_LE,gumbel_pred_U2_LE)
  
  hist_true <- hist2d(U1_LE, U2_LE, nbins = c(10,10), show = FALSE)
  
  hist_gauss <- hist2d(gauss_pred_U1_LE, gauss_pred_U2_LE, nbins = c(10,10), show = FALSE)
  hist_t <- hist2d(t_pred_U1_LE, t_pred_U2_LE, nbins = c(10,10), show = FALSE)
  hist_clayton <- hist2d(clayton_pred_U1_LE, clayton_pred_U2_LE, nbins = c(10,10), show = FALSE)
  hist_gumbel <- hist2d(gumbel_pred_U1_LE, gumbel_pred_U2_LE, nbins = c(10,10), show = FALSE)
  hist_frank <- hist2d(frank_pred_U1_LE, frank_pred_U2_LE, nbins = c(10,10), show = FALSE)
  
  par(mar = c(1,1,1,1), mfrow = c(2,3))
  
  hist3D(
    x = hist_true$x,
    y = hist_true$y,
    z = hist_true$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "Observed copula"
  )
  
  hist3D(
    x = hist_gauss$x,
    y = hist_gauss$y,
    z = hist_gauss$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "Gauss"
  )
  
  hist3D(
    x = hist_t$x,
    y = hist_t$y,
    z = hist_t$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "student t"
  )
  
  hist3D(
    x = hist_clayton$x,
    y = hist_clayton$y,
    z = hist_clayton$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "Clayton"
  )
  
  hist3D(
    x = hist_frank$x,
    y = hist_frank$y,
    z = hist_frank$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "Frank"
  )
  
  hist3D(
    x = hist_gumbel$x,
    y = hist_gumbel$y,
    z = hist_gumbel$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "gumbel"
  )
  
  ################################################################################
  # helper for test statistics
  
  library(cramer)
  
  cramer.test(cbind(U1_LE,U2_LE), cbind(U1_LE,U2_LE))
  cramer.test(cbind(U1_LE,U2_LE), cbind(t_pred_U1_LE, t_pred_U2_LE), replicates = 10000, sim = "permutation")
  cramer.test(cbind(U1_LE,U2_LE), cbind(gauss_pred_U1_LE, gauss_pred_U2_LE), replicates = 10000)
  cramer.test(cbind(U1_LE,U2_LE), cbind(frank_pred_U1_LE, frank_pred_U2_LE), replicates = 10000)
  cramer.test(cbind(U1_LE,U2_LE), cbind(clayton_pred_U1_LE, clayton_pred_U2_LE), replicates = 10000)
  cramer.test(cbind(U1_LE,U2_LE), cbind(gumbel_pred_U1_LE, gumbel_pred_U2_LE), replicates = 10000)
  
  library(fasano.franceschini.test)
  
  fasano.franceschini.test(cbind(U1_LE,U2_LE), cbind(U1_LE,U2_LE))
  fasano.franceschini.test(cbind(U1_LE,U2_LE), cbind(t_pred_U1_LE, t_pred_U2_LE))
  fasano.franceschini.test(cbind(U1_LE,U2_LE), cbind(gauss_pred_U1_LE, gauss_pred_U2_LE))
  fasano.franceschini.test(cbind(U1_LE,U2_LE), cbind(frank_pred_U1_LE, frank_pred_U2_LE))
  fasano.franceschini.test(cbind(U1_LE,U2_LE), cbind(clayton_pred_U1_LE, clayton_pred_U2_LE))
  fasano.franceschini.test(cbind(U1_LE,U2_LE), cbind(gumbel_pred_U1_LE, gumbel_pred_U2_LE))
  
  
}