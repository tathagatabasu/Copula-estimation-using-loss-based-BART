library(dplyr)

cia_wf_data <- read.csv("countries.csv")

colnames_for_anlysis <- c("Country","People.and.Society..Life.expectancy.at.birth...male",
                          "People.and.Society..Life.expectancy.at.birth...female",
                          "Economy..Real.GDP..purchasing.power.parity.")

cia_wf_data_le_vs_gdp <- cia_wf_data %>% dplyr::select(all_of(colnames_for_anlysis))

index_2023 <- (grepl("(2023 est.)", cia_wf_data_le_vs_gdp$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))


cia_wf_data_2023 <- cia_wf_data_le_vs_gdp[index_2023,]

index_bil <- (grepl("billion (2023 est.)", cia_wf_data_2023$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))
index_tril <- (grepl("trillion (2023 est.)", cia_wf_data_2023$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))
index_mil <- (grepl("million (2023 est.)", cia_wf_data_2023$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))


colnames(cia_wf_data_2023) <- c("Country",
                              "Life_expectancy_M",
                              "Life_expectancy_F",
                              "GDP_PPP")

library(readr)

cia_wf_data_2023 = cia_wf_data_2023 %>%
  mutate(across(-c(Country),.fns = parse_number))

# commented as repating this will be problematic

# cia_wf_data_2023$GDP_PPP[index_bil] = 1000 * cia_wf_data_2023$GDP_PPP[index_bil]
# cia_wf_data_2023$GDP_PPP[index_tril] = 1000000 * cia_wf_data_2023$GDP_PPP[index_tril]

# cia_wf_data_2023 <- na.omit(cia_wf_data_2023)

plot(log(cia_wf_data_2023$GDP_PPP),cia_wf_data_2023$Life_expectancy_M)
plot(log(cia_wf_data_2023$GDP_PPP),cia_wf_data_2023$Life_expectancy_F)

U1 = ecdf(cia_wf_data_2023$Life_expectancy_F)(cia_wf_data_2023$Life_expectancy_F)
U2 = ecdf(cia_wf_data_2023$Life_expectancy_M)(cia_wf_data_2023$Life_expectancy_M)

plot(U1,U2)
plot(cia_wf_data_2023$Life_expectancy_F,cia_wf_data_2023$Life_expectancy_M)

GDP <- as.data.frame((log(cia_wf_data_2023$GDP_PPP) - min(log(cia_wf_data_2023$GDP_PPP)))/(max(log(cia_wf_data_2023$GDP_PPP)) - min(log(cia_wf_data_2023$GDP_PPP))))
GDP <- as.matrix(GDP)
rownames(GDP) <- 1:nrow(GDP)

n.chain_par <- 1
n.iter_par <- 6000
n.born.out.par <- 1000
n.thin <- 1
incl.split_par <- TRUE
cont.unif_par <- TRUE
moves.prob_par <- c(0.4, 0.4, 0.1, 0.1)
##########################################################

gauss_GDP_tree_1 <- MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                 n.tree = 5,
                                 X = GDP,
                                 U1 = U1,
                                 U2 = U2,
                                 prior_list = lb.prior.def, 
                                 moves.prob = moves.prob_par, 
                                 starting.tree = NULL,
                                 cont.unif = cont.unif_par,
                                 include.split = incl.split_par,
                                 prop_mu = 0, prop_sigma = .2,
                                 theta_param_1 = 0, theta_param_2 = .3,
                                 var_param_1 = 1, var_param_2 = 2,
                                 prior_type = "N",
                                 cop_type = "gauss")

t_GDP_tree_1 <- MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                n.tree = 5,
                                X = GDP,
                                U1 = U1,
                                U2 = U2,
                                prior_list = lb.prior.def, 
                                moves.prob = moves.prob_par, 
                                starting.tree = NULL,
                                cont.unif = cont.unif_par,
                                include.split = incl.split_par,
                                prop_mu = 0, prop_sigma = .2,
                                theta_param_1 = 0, theta_param_2 = .3,
                                var_param_1 = 1, var_param_2 = 2,
                                prior_type = "N",
                                cop_type = "t")


model_twin <- t_GDP_tree_1

list_pred_lb <- lapply(1:length(model_twin$trees), \(idx) BART_calculate_pred(model_twin$trees[[idx]], GDP))

pred_val = do.call(rbind,list_pred_lb)

n.thin <- 1
n.iter_par <- 6000
n.born.out.par <- 1000

pred_val_vec = as.vector(pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])

pred_obs = rep(GDP, each = (n.chain_par * (n.iter_par - n.born.out.par)))

# like

like_val <- apply(pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1, U2))

like_df <-data.frame("nn" = like_val)
like_df$idx <- 1:(n.chain_par*n.iter_par)

pl_like <- ggplot(like_df, aes(idx, nn)) + 
  geom_line() + 
  ylab('log-likelihood') +
  theme_classic() + 
  theme(panel.grid.major = element_line())

pl_like


# theta_true = rep(param_gauss(get(paste0("tau_true_pred_",test_case))), each = (n.chain_par * (n.iter_par - n.born.out.par)))

pred_cond <- data.frame("obs" = pred_obs)
pred_cond$obs = pred_obs
# pred_cond$theta_true = theta_true
pred_cond$y = link_gauss(pred_val_vec)
pred_cond$tau = BiCopPar2Tau(1, pred_cond$y)

pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])

pred_cond_mod = pred_cond_thin %>%
  group_by(obs) %>%
  summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 

ggplot(pred_cond_mod) +
  geom_line(aes(obs, theta_mean)) +
  # geom_line(aes(obs, theta_true), col = 2) +
  geom_line(aes(obs, theta_q975), col = 3) +
  geom_line(aes(obs, theta_q025), col = 3) +
  # facet_wrap(facets = ~panel.name, ncol = 2) +
  xlab('X') +
  ylab('estimated rho') +
  theme_classic()

pred_cond_mod_tau = pred_cond_thin %>%
  group_by(obs) %>%
  summarise(tau_mean = mean(tau), tau_q975 = quantile(tau, .975), tau_q025 = quantile(tau, .025)) 

ggplot(pred_cond_mod_tau) +
  geom_line(aes(obs, tau_mean)) +
  # geom_line(aes(obs, theta_true), col = 2) +
  geom_line(aes(obs, tau_q975), col = 3) +
  geom_line(aes(obs, tau_q025), col = 3) +
  # facet_wrap(facets = ~panel.name, ncol = 2) +
  xlab('X') +
  ylab('estimated tau') +
  theme_classic()

GDP_gauss_pred <- BiCopSim(N = nrow(GDP), family = 1, par = link_gauss(colMeans(pred_val)))

pred_U1 = GDP_gauss_pred[,1]
pred_U2 = GDP_gauss_pred[,2]

plot(U1,U2)
plot(pred_U1,pred_U2)
