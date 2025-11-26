# packages
library(dplyr)
library(readr)
library(VineCopula)
library(ggplot2)
library(plot3D)
library(gplots)
library(patchwork)

set.seed(1e3)

# dataset
cia_wf_data <- read.csv("countries.csv")

cia_wf_data <- cia_wf_data %>% dplyr::select(all_of(c("Country","People.and.Society..Life.expectancy.at.birth...male",
                                                      "People.and.Society..Life.expectancy.at.birth...female",
                                                      "People.and.Society..Literacy...male",
                                                      "People.and.Society..Literacy...female",
                                                      "Economy..Real.GDP.per.capita")))
# "Economy..Real.GDP.per.capita")))

colnames(cia_wf_data) <- c("Country",
                           "Life_expectancy_M",
                           "Life_expectancy_F",
                           "Liter_M",
                           "Liter_F",
                           "GDP_PC")

cia_wf_data = cia_wf_data %>%
  mutate(across(-c(Country),.fns = parse_number))

cia_wf_data <- na.omit(cia_wf_data)

plot(log(cia_wf_data$GDP_PC),cia_wf_data$Life_expectancy_M)
plot(log(cia_wf_data$GDP_PC),cia_wf_data$Life_expectancy_F)

plot(log(cia_wf_data$GDP_PC),cia_wf_data$Liter_M)
plot(log(cia_wf_data$GDP_PC),cia_wf_data$Liter_F)

U1_LE = pobs(cia_wf_data$Life_expectancy_F)
U2_LE = pobs(cia_wf_data$Life_expectancy_M)

U1_LT = pobs(cia_wf_data$Liter_F)
U2_LT = pobs(cia_wf_data$Liter_M)

plot(U1_LE,U2_LE)
plot(cia_wf_data$Life_expectancy_F,cia_wf_data$Life_expectancy_M)

plot(U1_LT,U2_LT)
plot(cia_wf_data$Liter_F,cia_wf_data$Liter_M)

GDP <- as.data.frame((log(cia_wf_data$GDP_PC) - min(log(cia_wf_data$GDP_PC)))/(max(log(cia_wf_data$GDP_PC)) - min(log(cia_wf_data$GDP_PC))))
GDP <- as.matrix(GDP)
rownames(GDP) <- 1:nrow(GDP)

n.chain_par <- 5
n.iter_par <- 25000
n.born.out.par <- 500
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

################################################################################

################################################################################
load(paste0("gauss_GDP_LE_tree_",n.tree, ".Rdata"))

gauss_pred_val <- do.call(rbind,lapply(1:length(gauss_GDP_LE$trees), \(idx) BART_calculate_pred(gauss_GDP_LE$trees[[idx]], GDP)))

gauss_like_df <-data.frame("nn" = apply(gauss_pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1_LE, U2_LE)))
gauss_like_df$idx <- rep(1:n.iter_par, n.chain_par)
gauss_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)

gauss_pl_like <- gauss_like_df %>%
  ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  geom_line() +
  labs(
    x = "Iteration",
    y = "Log-likelihood"
  ) +
  guides(color = "none") +
  theme_minimal()

save(gauss_pl_like, file = "plot_gauss_LE.Rdata")
save(gauss_pred_val, file = "dat_gauss_LE.Rdata")

rm(gauss_pred_val, gauss_like_df, gauss_pl_like, gauss_GDP_LE)

gc()
gc()

load(paste0("gauss_GDP_LE_adapt_tree_",n.tree, ".Rdata"))

gauss_pred_val <- do.call(rbind,lapply(1:length(gauss_GDP_LE_adapt$trees), \(idx) BART_calculate_pred(gauss_GDP_LE_adapt$trees[[idx]], GDP)))

gauss_like_df <-data.frame("nn" = apply(gauss_pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1_LE, U2_LE)))
gauss_like_df$idx <- rep(1:n.iter_par, n.chain_par)
gauss_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)

gauss_pl_like <- gauss_like_df %>%
  ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  geom_line() +
  labs(
    x = "Iteration",
    y = "Log-likelihood"
  ) +
  guides(color = "none") +
  theme_minimal()

save(gauss_pl_like, file = "plot_gauss_LE_adapt.Rdata")
save(gauss_pred_val, file = "dat_gauss_LE_adapt.Rdata")

rm(gauss_pred_val, gauss_like_df, gauss_pl_like, gauss_GDP_LE_adapt)

gc()
gc()

load(paste0("gauss_GDP_LT_tree_",n.tree, ".Rdata"))

gauss_pred_val <- do.call(rbind,lapply(1:length(gauss_GDP_LT$trees), \(idx) BART_calculate_pred(gauss_GDP_LT$trees[[idx]], GDP)))

rm(gauss_GDP_LT)
gc()

gauss_like_df <-data.frame("nn" = apply(gauss_pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1_LT, U2_LT)))
gauss_like_df$idx <- rep(1:n.iter_par, n.chain_par)
gauss_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)

gauss_pl_like <- gauss_like_df %>%
  ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  geom_line() +
  labs(
    x = "Iteration",
    y = "Log-likelihood"
  ) +
  guides(color = "none") +
  theme_minimal()

save(gauss_pl_like, file = "plot_gauss_LT.Rdata")
save(gauss_pred_val, file = "dat_gauss_LT.Rdata")

rm(gauss_pred_val, gauss_like_df, gauss_pl_like, gauss_GDP_LT)

gc()
gc()

load(paste0("gauss_GDP_LT_adapt_tree_",n.tree, ".Rdata"))

gauss_pred_val <- do.call(rbind,lapply(1:length(gauss_GDP_LT_adapt$trees), \(idx) BART_calculate_pred(gauss_GDP_LT_adapt$trees[[idx]], GDP)))

rm(gauss_GDP_LT_adapt)
gc()

gauss_like_df <-data.frame("nn" = apply(gauss_pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1_LT, U2_LT)))
gauss_like_df$idx <- rep(1:n.iter_par, n.chain_par)
gauss_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)

gauss_pl_like <- gauss_like_df %>%
  ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  geom_line() +
  labs(
    x = "Iteration",
    y = "Log-likelihood"
  ) +
  guides(color = "none") +
  theme_minimal()

save(gauss_pl_like, file = "plot_gauss_LT_adapt.Rdata")
save(gauss_pred_val, file = "dat_gauss_LT_adapt.Rdata")

rm(gauss_pred_val, gauss_like_df, gauss_pl_like, gauss_GDP_LT_adapt)

gc()
gc()

load(paste0("t_GDP_LE_tree_",n.tree, ".Rdata"))

t_pred_val <- do.call(rbind,lapply(1:length(t_GDP_LE$trees), \(idx) BART_calculate_pred(t_GDP_LE$trees[[idx]], GDP)))

rm(t_GDP_LE)
gc()

t_like_df <-data.frame("nn" = apply(t_pred_val, 1, function(x)loglik_t(link_t(x), U1_LE, U2_LE)))
t_like_df$idx <- rep(1:n.iter_par, n.chain_par)
t_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)

t_pl_like <- t_like_df %>%
  ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  geom_line() +
  labs(
    x = "Iteration",
    y = "Log-likelihood"
  ) +
  guides(color = "none") +
  theme_minimal()

save(t_pl_like, file = "plot_t_LE.Rdata")
save(t_pred_val, file = "dat_t_LE.Rdata")

rm(t_pred_val, t_like_df, t_pl_like, t_GDP_LE)

gc()
gc()

load(paste0("t_GDP_LE_adapt_tree_",n.tree, ".Rdata"))

t_pred_val <- do.call(rbind,lapply(1:length(t_GDP_LE_adapt$trees), \(idx) BART_calculate_pred(t_GDP_LE_adapt$trees[[idx]], GDP)))

rm(t_GDP_LE_adapt)
gc()

t_like_df <-data.frame("nn" = apply(t_pred_val, 1, function(x)loglik_t(link_t(x), U1_LE, U2_LE)))
t_like_df$idx <- rep(1:n.iter_par, n.chain_par)
t_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)

t_pl_like <- t_like_df %>%
  ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  geom_line() +
  labs(
    x = "Iteration",
    y = "Log-likelihood"
  ) +
  guides(color = "none") +
  theme_minimal()

save(t_pl_like, file = "plot_t_LE_adapt.Rdata")
save(t_pred_val, file = "dat_t_LE_adapt.Rdata")

rm(t_pred_val, t_like_df, t_pl_like, t_GDP_LE_adapt)

gc()
gc()

load(paste0("t_GDP_LT_tree_",n.tree, ".Rdata"))

t_pred_val <- do.call(rbind,lapply(1:length(t_GDP_LT$trees), \(idx) BART_calculate_pred(t_GDP_LT$trees[[idx]], GDP)))

rm(t_GDP_LT)
gc()

t_like_df <-data.frame("nn" = apply(t_pred_val, 1, function(x)loglik_t(link_t(x), U1_LT, U2_LT)))
t_like_df$idx <- rep(1:n.iter_par, n.chain_par)
t_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)

t_pl_like <- t_like_df %>%
  ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  geom_line() +
  labs(
    x = "Iteration",
    y = "Log-likelihood"
  ) +
  guides(color = "none") +
  theme_minimal()

save(t_pl_like, file = "plot_t_LT.Rdata")
save(t_pred_val, file = "dat_t_LT.Rdata")

rm(t_pred_val, t_like_df, t_pl_like, t_GDP_LT)

gc()
gc()

load(paste0("t_GDP_LT_adapt_tree_",n.tree, ".Rdata"))

t_pred_val <- do.call(rbind,lapply(1:length(t_GDP_LT_adapt$trees), \(idx) BART_calculate_pred(t_GDP_LT_adapt$trees[[idx]], GDP)))

rm(t_GDP_LT_adapt)
gc()

t_like_df <-data.frame("nn" = apply(t_pred_val, 1, function(x)loglik_t(link_t(x), U1_LT, U2_LT)))
t_like_df$idx <- rep(1:n.iter_par, n.chain_par)
t_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)

t_pl_like <- t_like_df %>%
  ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  geom_line() +
  labs(
    x = "Iteration",
    y = "Log-likelihood"
  ) +
  guides(color = "none") +
  theme_minimal()

save(t_pl_like, file = "plot_t_LT_adapt.Rdata")
save(t_pred_val, file = "dat_t_LT_adapt.Rdata")

rm(t_pred_val, t_like_df, t_pl_like, t_GDP_LT_adapt)

gc()
gc()

################################################################################

################################################################################

if(F){
  
  if(n.tree == 5){
    load("real_case_5_trees/dat_gauss_LE.Rdata")
    load("real_case_5_trees/dat_t_LE.Rdata")
  }
  
  if(n.tree == 10){
    load("10_trees/dat_gauss_LE.Rdata")
    load("10_trees/dat_t_LE.Rdata")
  }
  
  # data plotting
  pred_cond <- data.frame("obs" = rep(GDP, each = (n.chain_par * (n.iter_par - n.born.out.par))))
  pred_cond$gauss_y = link_gauss(as.vector(gauss_pred_val))
  pred_cond$gauss_tau = BiCopPar2Tau(1, pred_cond$gauss_y)
  pred_cond$t_y = link_t(as.vector(t_pred_val))
  pred_cond$t_tau = BiCopPar2Tau(2, pred_cond$t_y)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),nrow(GDP))
  
  if(n.tree == 5){
    load("real_case_5_trees/dat_gauss_LE_adapt.Rdata")
    load("real_case_5_trees/dat_t_LE_adapt.Rdata")
  }
  
  if(n.tree == 10){
    load("10_trees/dat_gauss_LE_adapt.Rdata")
    load("10_trees/dat_t_LE_adapt.Rdata")
  }
  
  pred_cond_adapt <- data.frame("obs" = rep(GDP, each = (n.chain_par * (n.iter_par - n.born.out.par))))
  pred_cond_adapt$gauss_y = link_gauss(as.vector(gauss_pred_val))
  pred_cond_adapt$gauss_tau = BiCopPar2Tau(1, pred_cond_adapt$gauss_y)
  pred_cond_adapt$t_y = link_t(as.vector(t_pred_val))
  pred_cond_adapt$t_tau = BiCopPar2Tau(2, pred_cond_adapt$t_y)
  pred_cond_adapt$idx = rep(rep(1:n.iter_par, n.chain_par),nrow(GDP))
  
  pred_cond_mod = pred_cond %>%
    filter(idx > 5000) %>%
    group_by(obs) %>%
    summarise(gauss_theta_mean = mean(gauss_y), gauss_theta_q975 = quantile(gauss_y, .975), gauss_theta_q025 = quantile(gauss_y, .025),
              t_theta_mean = mean(t_y), t_theta_q975 = quantile(t_y, .975), t_theta_q025 = quantile(t_y, .025))
  
  pred_cond_mod_adapt = pred_cond_adapt %>%
    filter(idx > 5000) %>%
    group_by(obs) %>%
    summarise(gauss_theta_mean = mean(gauss_y), gauss_theta_q975 = quantile(gauss_y, .975), gauss_theta_q025 = quantile(gauss_y, .025),
              t_theta_mean = mean(t_y), t_theta_q975 = quantile(t_y, .975), t_theta_q025 = quantile(t_y, .025))
  
  hist_true <- hist2d(U1_LE, U2_LE, nbins = c(10,10), show = FALSE)
  
  set.seed(1e3)
  
  GDP_gauss_pred <- BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod$gauss_theta_mean)
  
  gauss_pred_U1_LE = GDP_gauss_pred[,1]
  gauss_pred_U2_LE = GDP_gauss_pred[,2]
  
  hist_gauss <- hist2d(gauss_pred_U1_LE, gauss_pred_U2_LE, nbins = c(10,10), show = FALSE)
  
  GDP_gauss_pred_adapt <- BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod_adapt$gauss_theta_mean)
  
  gauss_pred_U1_LE_adapt = GDP_gauss_pred_adapt[,1]
  gauss_pred_U2_LE_adapt = GDP_gauss_pred_adapt[,2]
  
  hist_gauss_adapt <- hist2d(gauss_pred_U1_LE_adapt, gauss_pred_U2_LE_adapt, nbins = c(10,10), show = FALSE)
  
  gc()
  
  GDP_t_pred <- BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod$t_theta_mean, par2 = 3)
  
  t_pred_U1_LE = GDP_t_pred[,1]
  t_pred_U2_LE = GDP_t_pred[,2]
  
  hist_t <- hist2d(t_pred_U1_LE, t_pred_U2_LE, nbins = c(10,10), show = FALSE)
  
  GDP_t_pred_adapt <- BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod_adapt$t_theta_mean, par2 = 3)
  
  t_pred_U1_LE_adapt = GDP_t_pred_adapt[,1]
  t_pred_U2_LE_adapt = GDP_t_pred_adapt[,2]
  
  hist_t_adapt <- hist2d(t_pred_U1_LE_adapt, t_pred_U2_LE_adapt, nbins = c(10,10), show = FALSE)
  
  gc()
  
  pred_cond_mod_tau = pred_cond %>%
    filter(idx > 5000) %>%
    group_by(obs) %>%
    summarise(gauss_tau_mean = mean(gauss_tau), gauss_tau_q975 = quantile(gauss_tau, .975), gauss_tau_q025 = quantile(gauss_tau, .025),
              t_tau_mean = mean(t_tau), t_tau_q975 = quantile(t_tau, .975), t_tau_q025 = quantile(t_tau, .025))
  
  pred_cond_mod_tau_adapt = pred_cond_adapt %>%
    filter(idx > 5000) %>%
    group_by(obs) %>%
    summarise(gauss_tau_mean = mean(gauss_tau), gauss_tau_q975 = quantile(gauss_tau, .975), gauss_tau_q025 = quantile(gauss_tau, .025),
              t_tau_mean = mean(t_tau), t_tau_q975 = quantile(t_tau, .975), t_tau_q025 = quantile(t_tau, .025))
  
  pred_cond_burn = pred_cond %>%
    filter(idx > 5000)
  
  pred_cond_adapt_burn = pred_cond_adapt %>%
    filter(idx > 5000)
  
  pl_tau_est <- ggplot(pred_cond_mod_tau) +
    geom_line(aes(obs, gauss_tau_mean, col = "gauss")) +
    geom_line(aes(obs, gauss_tau_q975, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, gauss_tau_q025, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, t_tau_mean, col = "t")) +
    geom_line(aes(obs, t_tau_q975, col = "t"),linetype="dotted") +
    geom_line(aes(obs, t_tau_q025, col = "t"),linetype="dotted") +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pl_tau_est
  
  pl_tau_est_adapt <- ggplot(pred_cond_mod_tau_adapt) +
    geom_line(aes(obs, gauss_tau_mean, col = "gauss")) +
    geom_line(aes(obs, gauss_tau_q975, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, gauss_tau_q025, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, t_tau_mean, col = "t")) +
    geom_line(aes(obs, t_tau_q975, col = "t"),linetype="dotted") +
    geom_line(aes(obs, t_tau_q025, col = "t"),linetype="dotted") +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pl_tau_est_adapt
  
  save(pred_cond_mod_tau_adapt, pred_cond_mod_tau, pl_tau_est, pl_tau_est_adapt, file = "pred_tau_LE_adapt.Rdata")
  
  cor(U1_LE,U2_LE, method = "kendall")
  
  par(mar = c(5,5,2,1), mfrow = c(2,3))
  
  plot(U1_LE,U2_LE, xlab = "F Life Exp", ylab = "M Life Exp", main = "Observed data")
  plot(gauss_pred_U1_LE_adapt,gauss_pred_U2_LE_adapt, xlab = "F Life Exp", ylab = "M Life Exp", main = "Predicted (Gaussian)")
  plot(t_pred_U1_LE_adapt,t_pred_U2_LE_adapt, xlab = "F Life Exp", ylab = "M Life Exp", main = "Predicted (Student-t)")
  
  hist3D(
    x = hist_true$x,
    y = hist_true$y,
    z = hist_true$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss_adapt$counts, hist_t_adapt$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Life Exp", ylab = "M Life Exp", zlab = ""
  )
  
  hist3D(
    x = hist_gauss_adapt$x,
    y = hist_gauss_adapt$y,
    z = hist_gauss_adapt$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss_adapt$counts, hist_t_adapt$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Life Exp", ylab = "M Life Exp", zlab = ""
  )
  
  hist3D(
    x = hist_t_adapt$x,
    y = hist_t_adapt$y,
    z = hist_t_adapt$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss_adapt$counts, hist_t_adapt$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Life Exp", ylab = "M Life Exp", zlab = ""
  )
  
  # 8 5.5
  
  plot(U1_LE,U2_LE, xlab = "F Life Exp", ylab = "M Life Exp", main = "Observed data")
  plot(gauss_pred_U1_LE,gauss_pred_U2_LE, xlab = "F Life Exp", ylab = "M Life Exp", main = "Predicted (Gaussian)")
  plot(t_pred_U1_LE,t_pred_U2_LE, xlab = "F Life Exp", ylab = "M Life Exp", main = "Predicted (Student-t)")
  hist3D(
    x = hist_true$x,
    y = hist_true$y,
    z = hist_true$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss$counts, hist_t$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Life Exp", ylab = "M Life Exp", zlab = ""
  )
  
  hist3D(
    x = hist_gauss$x,
    y = hist_gauss$y,
    z = hist_gauss$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss$counts, hist_t$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Life Exp", ylab = "M Life Exp", zlab = ""
  )
  
  hist3D(
    x = hist_t$x,
    y = hist_t$y,
    z = hist_t$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss$counts, hist_t$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Life Exp", ylab = "M Life Exp", zlab = ""
  )
    # helper for test statistics
    
    library(cramer)
    
    cram_gauss <- sapply(1:100, function(i){set.seed(i); cramer.test(cbind(U1_LE,U2_LE), BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod$gauss_theta_mean), replicates = 1000, sim = "permutation")$p.value})
    cram_t <- sapply(1:100, function(i){set.seed(i);cramer.test(cbind(U1_LE,U2_LE), BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod$t_theta_mean, par2 = 3), replicates = 1000)$p.value})
    
    cram_gauss_adapt <- sapply(1:100, function(i){set.seed(i); cramer.test(cbind(U1_LE,U2_LE), BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod_adapt$gauss_theta_mean), replicates = 1000, sim = "permutation")$p.value})
    cram_t_adapt <- sapply(1:100, function(i){set.seed(i);cramer.test(cbind(U1_LE,U2_LE), BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod_adapt$t_theta_mean, par2 = 3), replicates = 1000)$p.value})
    
    library(fasano.franceschini.test)
    
    ff_gauss <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(cbind(U1_LE,U2_LE), BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod$gauss_theta_mean), nPermute = 1000, verbose = F)$p.value})
    ff_t <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(cbind(U1_LE,U2_LE), BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod$t_theta_mean, par2 = 3), nPermute = 1000, verbose = F)$p.value})
    
    ff_gauss_adapt <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(cbind(U1_LE,U2_LE), BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod_adapt$gauss_theta_mean), nPermute = 1000, verbose = F)$p.value})
    ff_t_adapt <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(cbind(U1_LE,U2_LE), BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod_adapt$t_theta_mean, par2 = 3), nPermute = 1000, verbose = F)$p.value})
    
    p_val_summ_adapt <- rbind(c(mean(cram_gauss), median(cram_gauss), sd(cram_gauss), mean(ff_gauss), median(ff_gauss), sd(ff_gauss)),
                              c(mean(cram_t), median(cram_t), sd(cram_t), mean(ff_t), median(ff_t), sd(ff_t)),
                              c(mean(cram_gauss_adapt), median(cram_gauss_adapt), sd(cram_gauss_adapt), mean(ff_gauss_adapt), median(ff_gauss_adapt), sd(ff_gauss_adapt)),
                              c(mean(cram_t_adapt), median(cram_t_adapt), sd(cram_t_adapt), mean(ff_t_adapt), median(ff_t_adapt), sd(ff_t_adapt)))
    
    library(xtable)
    xtable(p_val_summ_adapt)
    
    pl_tau_est + ylim(0,1) + geom_hline(yintercept = cor(U1_LE,U2_LE, method = "kendall"), linetype = "dotted", linewidth = 0.2) + labs(title="Without adaption") + xlab("Scaled GDP") + pl_tau_est_adapt + ylim(0,1) + geom_hline(yintercept = cor(U1_LE,U2_LE, method = "kendall"), linetype = "dotted", linewidth = 0.2) + labs(title="With adaption") + xlab("Scaled GDP")
}

if(F){
  
  gauss_woa_10 <- get(load("10_trees/plot_gauss_LE.Rdata"))
  gauss_wa_10 <- get(load("10_trees/plot_gauss_LE_adapt.Rdata"))
  
  t_woa_10 <- get(load("10_trees/plot_t_LE.Rdata"))
  t_wa_10 <- get(load("10_trees/plot_t_LE_adapt.Rdata"))
  
  gauss_woa_5 <- get(load("real_case_5_trees/plot_gauss_LE.Rdata"))
  gauss_wa_5 <- get(load("real_case_5_trees/plot_gauss_LE_adapt.Rdata"))
  
  t_woa_5 <- get(load("real_case_5_trees/plot_t_LE.Rdata"))
  t_wa_5 <- get(load("real_case_5_trees/plot_t_LE_adapt.Rdata"))
  
  (gauss_woa_10 + ylim(200,280) + labs(title="Gaussian (without adaption)") + gauss_wa_10 + ylim(200,280) + labs(title="Gaussian (with adaption)")) / 
    (t_woa_10 + ylim(200,280) + labs(title="Student-t (without adaption)") + t_wa_10 + ylim(200,280) + labs(title="Student-t (with adaption)"))
  
  (gauss_woa_5 + ylim(100,260) + labs(title="Gaussian (without adaption)") + gauss_wa_5 + ylim(100,260) + labs(title="Gaussian (with adaption)")) / 
    (t_woa_5 + ylim(120,260) + labs(title="Student-t (without adaption)") + t_wa_5 + ylim(120,260) + labs(title="Student-t (with adaption)")) #/
  
  # par(mar = c(5,5,5,1), mfrow = c(2,2))
  # acf(mcmc(pred_cond_burn$gauss_y), lag.max = 200, main = "Gaussian (without adaption)")
  # acf(mcmc(pred_cond_adapt_burn$gauss_y), lag.max = 200, main = "Gaussian (with adaption)")
  # acf(mcmc(pred_cond_burn$t_y), lag.max = 200, main = "Student-t (without adaption)")
  # acf(mcmc(pred_cond_adapt_burn$t_y), lag.max = 200, main = "Student-t (with adaption)")
}


################################################################################

################################################################################

if(F){
  
  if(n.tree == 5){
    load("real_case_5_trees/dat_gauss_LT.Rdata")
    load("real_case_5_trees/dat_t_LT.Rdata")
  }
  
  if(n.tree == 5){
    load("real_case_5_trees/dat_gauss_LT_adapt.Rdata")
    load("real_case_5_trees/dat_t_LT_adapt.Rdata")
  }
  
  if(n.tree == 10){
    load("10_trees/dat_gauss_LT.Rdata")
    load("10_trees/dat_t_LT.Rdata")
  }
  
  if(n.tree == 10){
    load("10_trees/dat_gauss_LT_adapt.Rdata")
    load("10_trees/dat_t_LT_adapt.Rdata")
  }
  
  # data plotting
  pred_cond <- data.frame("obs" = rep(GDP, each = (n.chain_par * (n.iter_par - n.born.out.par))))
  pred_cond$gauss_y = link_gauss(as.vector(gauss_pred_val))
  pred_cond$gauss_tau = BiCopPar2Tau(1, pred_cond$gauss_y)
  pred_cond$t_y = link_t(as.vector(t_pred_val))
  pred_cond$t_tau = BiCopPar2Tau(2, pred_cond$t_y)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),nrow(GDP))
  
  pred_cond_adapt <- data.frame("obs" = rep(GDP, each = (n.chain_par * (n.iter_par - n.born.out.par))))
  pred_cond_adapt$gauss_y = link_gauss(as.vector(gauss_pred_val))
  pred_cond_adapt$gauss_tau = BiCopPar2Tau(1, pred_cond_adapt$gauss_y)
  pred_cond_adapt$t_y = link_t(as.vector(t_pred_val))
  pred_cond_adapt$t_tau = BiCopPar2Tau(2, pred_cond_adapt$t_y)
  pred_cond_adapt$idx = rep(rep(1:n.iter_par, n.chain_par),nrow(GDP))
  
  pred_cond_mod = pred_cond %>%
    filter(idx > 5000) %>%
    group_by(obs) %>%
    summarise(gauss_theta_mean = mean(gauss_y), gauss_theta_q975 = quantile(gauss_y, .975), gauss_theta_q025 = quantile(gauss_y, .025),
              t_theta_mean = mean(t_y), t_theta_q975 = quantile(t_y, .975), t_theta_q025 = quantile(t_y, .025))
  
  pred_cond_mod_adapt = pred_cond_adapt %>%
    filter(idx > 5000) %>%
    group_by(obs) %>%
    summarise(gauss_theta_mean = mean(gauss_y), gauss_theta_q975 = quantile(gauss_y, .975), gauss_theta_q025 = quantile(gauss_y, .025),
              t_theta_mean = mean(t_y), t_theta_q975 = quantile(t_y, .975), t_theta_q025 = quantile(t_y, .025))
  
  hist_true <- hist2d(U1_LT, U2_LT, nbins = c(10,10), show = FALSE)
  
  GDP_gauss_pred <- BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod$gauss_theta_mean)
  
  gauss_pred_U1_LT = GDP_gauss_pred[,1]
  gauss_pred_U2_LT = GDP_gauss_pred[,2]
  
  hist_gauss <- hist2d(gauss_pred_U1_LT, gauss_pred_U2_LT, nbins = c(10,10), show = FALSE)
  
  GDP_gauss_pred_adapt <- BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod_adapt$gauss_theta_mean)
  
  gauss_pred_U1_LT_adapt = GDP_gauss_pred_adapt[,1]
  gauss_pred_U2_LT_adapt = GDP_gauss_pred_adapt[,2]
  
  hist_gauss_adapt <- hist2d(gauss_pred_U1_LT_adapt, gauss_pred_U2_LT_adapt, nbins = c(10,10), show = FALSE)
  
  gc()
  
  GDP_t_pred <- BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod$t_theta_mean, par2 = 3)
  
  t_pred_U1_LT = GDP_t_pred[,1]
  t_pred_U2_LT = GDP_t_pred[,2]
  
  hist_t <- hist2d(t_pred_U1_LT, t_pred_U2_LT, nbins = c(10,10), show = FALSE)
  
  GDP_t_pred_adapt <- BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod_adapt$t_theta_mean, par2 = 3)
  
  t_pred_U1_LT_adapt = GDP_t_pred_adapt[,1]
  t_pred_U2_LT_adapt = GDP_t_pred_adapt[,2]
  
  hist_t_adapt <- hist2d(t_pred_U1_LT_adapt, t_pred_U2_LT_adapt, nbins = c(10,10), show = FALSE)
  
  gc()
  
  pred_cond_mod_tau = pred_cond %>%
    filter(idx > 5000) %>%
    group_by(obs) %>%
    summarise(gauss_tau_mean = mean(gauss_tau), gauss_tau_q975 = quantile(gauss_tau, .975), gauss_tau_q025 = quantile(gauss_tau, .025),
              t_tau_mean = mean(t_tau), t_tau_q975 = quantile(t_tau, .975), t_tau_q025 = quantile(t_tau, .025))
  
  pred_cond_mod_tau_adapt = pred_cond_adapt %>%
    filter(idx > 5000) %>%
    group_by(obs) %>%
    summarise(gauss_tau_mean = mean(gauss_tau), gauss_tau_q975 = quantile(gauss_tau, .975), gauss_tau_q025 = quantile(gauss_tau, .025),
              t_tau_mean = mean(t_tau), t_tau_q975 = quantile(t_tau, .975), t_tau_q025 = quantile(t_tau, .025))
  
  pred_cond_burn = pred_cond %>%
    filter(idx>5000)
  
  pred_cond_adapt_burn = pred_cond_adapt %>%
    filter(idx>5000)
  
  pl_tau_est <- ggplot(pred_cond_mod_tau) +
    geom_line(aes(obs, gauss_tau_mean, col = "gauss")) +
    geom_line(aes(obs, gauss_tau_q975, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, gauss_tau_q025, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, t_tau_mean, col = "t")) +
    geom_line(aes(obs, t_tau_q975, col = "t"),linetype="dotted") +
    geom_line(aes(obs, t_tau_q025, col = "t"),linetype="dotted") +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pl_tau_est
  
  pl_tau_est_adapt <- ggplot(pred_cond_mod_tau_adapt) +
    geom_line(aes(obs, gauss_tau_mean, col = "gauss")) +
    geom_line(aes(obs, gauss_tau_q975, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, gauss_tau_q025, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, t_tau_mean, col = "t")) +
    geom_line(aes(obs, t_tau_q975, col = "t"),linetype="dotted") +
    geom_line(aes(obs, t_tau_q025, col = "t"),linetype="dotted") +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pl_tau_est_adapt
  
  save(pred_cond_mod_tau_adapt, pred_cond_mod_tau, pl_tau_est, pl_tau_est_adapt, file = "pred_tau_LT_adapt.Rdata")
  
  cor(U1_LT,U2_LT, method = "kendall")
  
  par(mar = c(5,5,2,1), mfrow = c(2,3))
  
  plot(U1_LT,U2_LT, xlab = "F Literacy", ylab = "M Literacy", main = "Observed data")
  plot(gauss_pred_U1_LT_adapt,gauss_pred_U2_LT_adapt, xlab = "F Literacy", ylab = "M Literacy", main = "Simulated (Gaussian)")
  plot(t_pred_U1_LT_adapt,t_pred_U2_LT_adapt, xlab = "F Literacy", ylab = "M Literacy", main = "Simulated (Student-t)")
  
  hist3D(
    x = hist_true$x,
    y = hist_true$y,
    z = hist_true$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss_adapt$counts, hist_t_adapt$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Literacy", ylab = "M Literacy", zlab = ""
  )
  
  hist3D(
    x = hist_gauss_adapt$x,
    y = hist_gauss_adapt$y,
    z = hist_gauss_adapt$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss_adapt$counts, hist_t_adapt$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Literacy", ylab = "M Literacy", zlab = ""
  )
  
  hist3D(
    x = hist_t_adapt$x,
    y = hist_t_adapt$y,
    z = hist_t_adapt$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss_adapt$counts, hist_t_adapt$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Literacy", ylab = "M Literacy", zlab = ""
  )
  
  # 8 5
  
  plot(U1_LT,U2_LT, xlab = "F Literacy", ylab = "M Literacy", main = "Observed data")
  plot(gauss_pred_U1_LT,gauss_pred_U2_LT, xlab = "F Literacy", ylab = "M Literacy", main = "Simulated (Gaussian)")
  plot(t_pred_U1_LT,t_pred_U2_LT, xlab = "F Literacy", ylab = "M Literacy", main = "Simulated (Student-t)")
  
  hist3D(
    x = hist_true$x,
    y = hist_true$y,
    z = hist_true$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss$counts, hist_t$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Literacy", ylab = "M Literacy", zlab = ""
  )
  
  hist3D(
    x = hist_gauss$x,
    y = hist_gauss$y,
    z = hist_gauss$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss$counts, hist_t$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Literacy", ylab = "M Literacy", zlab = ""
  )
  
  hist3D(
    x = hist_t$x,
    y = hist_t$y,
    z = hist_t$counts,
    zlim = c(0, max(hist_true$counts, hist_gauss$counts, hist_t$counts)),
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.05, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "simple",
    xlab = "F Literacy", ylab = "M Literacy", zlab = ""
  )
  # helper for test statistics
  
  library(cramer)
  
  cram_gauss <- sapply(1:100, function(i){set.seed(i); cramer.test(cbind(U1_LT,U2_LT), BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod$gauss_theta_mean), replicates = 1000, sim = "permutation")$p.value})
  cram_t <- sapply(1:100, function(i){set.seed(i);cramer.test(cbind(U1_LT,U2_LT), BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod$t_theta_mean, par2 = 3), replicates = 1000)$p.value})
  
  cram_gauss_adapt <- sapply(1:100, function(i){set.seed(i); cramer.test(cbind(U1_LT,U2_LT), BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod_adapt$gauss_theta_mean), replicates = 1000, sim = "permutation")$p.value})
  cram_t_adapt <- sapply(1:100, function(i){set.seed(i);cramer.test(cbind(U1_LT,U2_LT), BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod_adapt$t_theta_mean, par2 = 3), replicates = 1000)$p.value})
  
  library(fasano.franceschini.test)
  
  ff_gauss <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(cbind(U1_LT,U2_LT), BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod$gauss_theta_mean), nPermute = 1000, verbose = F)$p.value})
  ff_t <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(cbind(U1_LT,U2_LT), BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod$t_theta_mean, par2 = 3), nPermute = 1000, verbose = F)$p.value})
  
  ff_gauss_adapt <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(cbind(U1_LT,U2_LT), BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod_adapt$gauss_theta_mean), nPermute = 1000, verbose = F)$p.value})
  ff_t_adapt <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(cbind(U1_LT,U2_LT), BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod_adapt$t_theta_mean, par2 = 3), nPermute = 1000, verbose = F)$p.value})
  
  p_val_summ_adapt <- rbind(c(mean(cram_gauss), median(cram_gauss), sd(cram_gauss), mean(ff_gauss), median(ff_gauss), sd(ff_gauss)),
                            c(mean(cram_t), median(cram_t), sd(cram_t), mean(ff_t), median(ff_t), sd(ff_t)),
                            c(mean(cram_gauss_adapt), median(cram_gauss_adapt), sd(cram_gauss_adapt), mean(ff_gauss_adapt), median(ff_gauss_adapt), sd(ff_gauss_adapt)),
                            c(mean(cram_t_adapt), median(cram_t_adapt), sd(cram_t_adapt), mean(ff_t_adapt), median(ff_t_adapt), sd(ff_t_adapt)))
  
  xtable(p_val_summ_adapt)
  
  pl_tau_est + ylim(-0.25,1) + geom_hline(yintercept = cor(U1_LT,U2_LT, method = "kendall"), linetype = "dotted", linewidth = 0.2) + labs(title="Without adaption") + xlab("Scaled GDP") + pl_tau_est_adapt + ylim(-0.25,1) + geom_hline(yintercept = cor(U1_LT,U2_LT, method = "kendall"), linetype = "dotted", linewidth = 0.2) + labs(title="With adaption") + xlab("Scaled GDP")
  
}

if(F){
  
  gauss_woa_10 <- get(load("10_trees/plot_gauss_LT.Rdata"))
  gauss_wa_10 <- get(load("10_trees/plot_gauss_LT_adapt.Rdata"))
  
  t_woa_10 <- get(load("10_trees/plot_t_LT.Rdata"))
  t_wa_10 <- get(load("10_trees/plot_t_LT_adapt.Rdata"))
  
  gauss_woa_5 <- get(load("real_case_5_trees/plot_gauss_LT.Rdata"))
  gauss_wa_5 <- get(load("real_case_5_trees/plot_gauss_LT_adapt.Rdata"))
  
  t_woa_5 <- get(load("real_case_5_trees/plot_t_LT.Rdata"))
  t_wa_5 <- get(load("real_case_5_trees/plot_t_LT_adapt.Rdata"))
  
  (gauss_woa_10 + ylim(310,360) + labs(title="Gaussian (without adaption)") + gauss_wa_10 + ylim(310,360) + labs(title="Gaussian (with adaption)")) / 
    (t_woa_10 + ylim(300,370) + labs(title="Student-t (without adaption)") + t_wa_10 + ylim(300,370) + labs(title="Student-t (with adaption)"))
  
  (gauss_woa_5 + ylim(310,350) + labs(title="Gaussian (without adaption)") + gauss_wa_5 + ylim(310,350) + labs(title="Gaussian (with adaption)")) / 
    (t_woa_5 + ylim(320,355) + labs(title="Student-t (without adaption)") + t_wa_5 + ylim(320,355) + labs(title="Student-t (with adaption)")) #/
  
  # par(mar = c(5,5,5,1), mfrow = c(2,2))
  # acf(mcmc(pred_cond_burn$gauss_y), lag.max = 200, main = "Gaussian (without adaption)")
  # acf(mcmc(pred_cond_adapt_burn$gauss_y), lag.max = 200, main = "Gaussian (with adaption)")
  # acf(mcmc(pred_cond_burn$t_y), lag.max = 200, main = "Student-t (without adaption)")
  # acf(mcmc(pred_cond_adapt_burn$t_y), lag.max = 200, main = "Student-t (with adaption)")
}
