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
library(patchwork)
library(dplyr)
library(readr)


set.seed(123)

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

plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Life_expectancy_M)
plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Life_expectancy_F)

plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Liter_M)
plot(log(cia_wf_data$GDP_PPP),cia_wf_data$Liter_F)

U1_LE = ecdf(cia_wf_data$Life_expectancy_F)(cia_wf_data$Life_expectancy_F)
U2_LE = ecdf(cia_wf_data$Life_expectancy_M)(cia_wf_data$Life_expectancy_M)

U1_LT = ecdf(cia_wf_data$Liter_F)(cia_wf_data$Liter_F)
U2_LT = ecdf(cia_wf_data$Liter_M)(cia_wf_data$Liter_M)

plot(U1_LE,U2_LE)
plot(cia_wf_data$Life_expectancy_F,cia_wf_data$Life_expectancy_M)

plot(U1_LT,U2_LT)
plot(cia_wf_data$Liter_F,cia_wf_data$Liter_M)

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

################################################################################

################################################################################

if(T){
  
  load("gauss_gdp_LE_tree_10.Rdata")
  
  gauss_list_pred_lb <- lapply(1:length(gauss_GDP_LE$trees), \(idx) BART_calculate_pred(gauss_GDP_LE$trees[[idx]], GDP))
  
  gauss_pred_val = do.call(rbind,gauss_list_pred_lb)
  
  gauss_pred_val_vec = as.vector(gauss_pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  gauss_like_val <- apply(gauss_pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1_LE, U2_LE))
  
  gauss_like_df <-data.frame("nn" = gauss_like_val)
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
  
  gauss_pl_like
  
  load("t_gdp_LE_tree_10.Rdata")
  
  t_list_pred_lb <- lapply(1:length(t_GDP_LE$trees), \(idx) BART_calculate_pred(t_GDP_LE$trees[[idx]], GDP))
  
  t_pred_val = do.call(rbind,t_list_pred_lb)
  
  t_pred_val_vec = as.vector(t_pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  t_like_val <- apply(t_pred_val, 1, function(x)loglik_t(link_t(x), U1_LE, U2_LE))
  
  t_like_df <-data.frame("nn" = t_like_val)
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
  
  t_pl_like
  
  
  load("clayton_gdp_LE_tree_10.Rdata")
  
  clayton_list_pred_lb <- lapply(1:length(clayton_GDP_LE$trees), \(idx) BART_calculate_pred(clayton_GDP_LE$trees[[idx]], GDP))
  
  clayton_pred_val = do.call(rbind,clayton_list_pred_lb)
  
  clayton_pred_val_vec = as.vector(clayton_pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  clayton_like_val <- apply(clayton_pred_val, 1, function(x)loglik_clayton(link_clayton(x), U1_LE, U2_LE))
  
  clayton_like_df <-data.frame("nn" = clayton_like_val)
  clayton_like_df$idx <- rep(1:n.iter_par, n.chain_par)
  clayton_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  clayton_pl_like <- clayton_like_df %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  clayton_pl_like
  
  load("gumbel_gdp_LE_tree_10.Rdata")
  
  gumbel_list_pred_lb <- lapply(1:length(gumbel_GDP_LE$trees), \(idx) BART_calculate_pred(gumbel_GDP_LE$trees[[idx]], GDP))
  
  gumbel_pred_val = do.call(rbind,gumbel_list_pred_lb)
  
  gumbel_pred_val_vec = as.vector(gumbel_pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  gumbel_like_val <- apply(gumbel_pred_val, 1, function(x)loglik_gumbel(link_gumbel(x), U1_LE, U2_LE))
  
  gumbel_like_df <-data.frame("nn" = gumbel_like_val)
  gumbel_like_df$idx <- rep(1:n.iter_par, n.chain_par)
  gumbel_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  gumbel_pl_like <- gumbel_like_df %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  gumbel_pl_like
  
  load("frank_gdp_LE_tree_10.Rdata")
  
  frank_list_pred_lb <- lapply(1:length(frank_GDP_LE$trees), \(idx) BART_calculate_pred(frank_GDP_LE$trees[[idx]], GDP))
  
  frank_pred_val = do.call(rbind,frank_list_pred_lb)
  
  frank_pred_val_vec = as.vector(frank_pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  frank_like_val <- apply(frank_pred_val, 1, function(x)loglik_frank(link_frank(x), U1_LE, U2_LE))
  
  frank_like_df <-data.frame("nn" = frank_like_val)
  frank_like_df$idx <- rep(1:n.iter_par, n.chain_par)
  frank_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  frank_pl_like <- frank_like_df %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  frank_pl_like
  
  
}

################################################################################

################################################################################

if(F){
  
  
  
  # data plotting
  pred_cond <- data.frame("obs" = rep(GDP, each = (n.chain_par * (n.iter_par - n.born.out.par))))
  pred_cond$gauss_y = link_gauss(gauss_pred_val_vec)
  pred_cond$gauss_tau = BiCopPar2Tau(1, pred_cond$gauss_y)
  pred_cond$t_y = link_t(t_pred_val_vec)
  pred_cond$t_tau = BiCopPar2Tau(2, pred_cond$t_y)
  pred_cond$clayton_y = link_clayton(clayton_pred_val_vec)
  pred_cond$clayton_tau = BiCopPar2Tau(3, pred_cond$clayton_y)
  pred_cond$gumbel_y = link_gumbel(gumbel_pred_val_vec)
  pred_cond$gumbel_tau = BiCopPar2Tau(4, pred_cond$gumbel_y)
  pred_cond$frank_y = link_frank(frank_pred_val_vec)
  pred_cond$frank_tau = BiCopPar2Tau(5, pred_cond$frank_y)
  
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(obs) %>%
    summarise(gauss_theta_mean = mean(gauss_y), gauss_theta_q975 = quantile(gauss_y, .975), gauss_theta_q025 = quantile(gauss_y, .025),
              frank_theta_mean = mean(frank_y), frank_theta_q975 = quantile(frank_y, .975), frank_theta_q025 = quantile(frank_y, .025),
              t_theta_mean = mean(t_y), t_theta_q975 = quantile(t_y, .975), t_theta_q025 = quantile(t_y, .025),
              clayton_theta_mean = mean(clayton_y), clayton_theta_q975 = quantile(clayton_y, .975), clayton_theta_q025 = quantile(clayton_y, .025),
              gumbel_theta_mean = mean(gumbel_y), gumbel_theta_q975 = quantile(gumbel_y, .975), gumbel_theta_q025 = quantile(gumbel_y, .025))
  
  GDP_gauss_pred <- BiCopSim(N = nrow(unique(GDP)), family = 1, par = pred_cond_mod$gauss_theta_mean)
  
  gauss_pred_U1_LT = GDP_gauss_pred[,1]
  gauss_pred_U2_LT = GDP_gauss_pred[,2]
  
  hist_gauss <- hist2d(gauss_pred_U1_LT, gauss_pred_U2_LT, nbins = c(10,10), show = FALSE)
  
  save(gauss_pred_U1_LT, gauss_pred_U2_LT, hist_gauss, gauss_like_val, gauss_pl_like, gauss_like_df, gauss_pred_val_vec, file = "gauss_gdp_LT_post.Rdata")
  rm(gauss_GDP_LT, gauss_pred_U1_LT, gauss_pred_U2_LT, hist_gauss, gauss_like_val, gauss_pl_like, gauss_like_df, GDP_gauss_pred, gauss_pred_val, gauss_list_pred_lb, gauss_pred_val_vec)
  gc()
  
  GDP_t_pred <- BiCopSim(N = nrow(unique(GDP)), family = 2, par = pred_cond_mod$t_theta_mean, par2 = 3)
  
  t_pred_U1_LT = GDP_t_pred[,1]
  t_pred_U2_LT = GDP_t_pred[,2]
  
  hist_t <- hist2d(t_pred_U1_LT, t_pred_U2_LT, nbins = c(10,10), show = FALSE)
  
  save(t_pred_U1_LT, t_pred_U2_LT, hist_t, t_like_val, t_pl_like, t_like_df, t_pred_val_vec, file = "t_gdp_LT_post.Rdata")
  rm(t_GDP_LT, t_pred_U1_LT, t_pred_U2_LT, hist_t, t_like_val, t_pl_like, t_like_df, GDP_t_pred, t_pred_val, t_list_pred_lb, t_pred_val_vec)
  gc()
  
  GDP_clayton_pred <- BiCopSim(N = nrow(unique(GDP)), family = 3, par = pred_cond_mod$clayton_theta_mean)
  
  clayton_pred_U1_LT = GDP_clayton_pred[,1]
  clayton_pred_U2_LT = GDP_clayton_pred[,2]
  
  hist_clayton <- hist2d(clayton_pred_U1_LT, clayton_pred_U2_LT, nbins = c(10,10), show = FALSE)
  
  save(clayton_pred_U1_LT, clayton_pred_U2_LT, hist_clayton, clayton_like_val, clayton_pl_like, clayton_like_df, clayton_pred_val_vec, file = "clayton_gdp_LT_post.Rdata")
  rm(clayton_GDP_LT, clayton_pred_U1_LT, clayton_pred_U2_LT, hist_clayton, clayton_like_val, clayton_pl_like, clayton_like_df, GDP_clayton_pred, clayton_pred_val, clayton_list_pred_lb, clayton_pred_val_vec)
  gc()
  
  GDP_gumbel_pred <- BiCopSim(N = nrow(unique(GDP)), family = 4, par = pred_cond_mod$gumbel_theta_mean)
  
  gumbel_pred_U1_LT = GDP_gumbel_pred[,1]
  gumbel_pred_U2_LT = GDP_gumbel_pred[,2]
  
  hist_gumbel <- hist2d(gumbel_pred_U1_LT, gumbel_pred_U2_LT, nbins = c(10,10), show = FALSE)
  
  save(gumbel_pred_U1_LT, gumbel_pred_U2_LT, hist_gumbel, gumbel_like_val, gumbel_pl_like, gumbel_like_df, gumbel_pred_val_vec, file = "gumbel_gdp_LT_post.Rdata")
  rm(gumbel_GDP_LT, gumbel_pred_U1_LT, gumbel_pred_U2_LT, hist_gumbel, gumbel_like_val, gumbel_pl_like, gumbel_like_df, GDP_gumbel_pred, gumbel_pred_val, gumbel_list_pred_lb, gumbel_pred_val_vec)
  gc()
  
  GDP_frank_pred <- BiCopSim(N = nrow(unique(GDP)), family = 5, par = pred_cond_mod$frank_theta_mean)
  
  frank_pred_U1_LT = GDP_frank_pred[,1]
  frank_pred_U2_LT = GDP_frank_pred[,2]
  
  hist_frank <- hist2d(frank_pred_U1_LT, frank_pred_U2_LT, nbins = c(10,10), show = FALSE)
  
  save(frank_pred_U1_LT, frank_pred_U2_LT, hist_frank, frank_like_val, frank_pl_like, frank_like_df, frank_pred_val_vec, file = "frank_gdp_LT_post.Rdata")
  rm(frank_GDP_LT, frank_pred_U1_LT, frank_pred_U2_LT, hist_frank, frank_like_val, frank_pl_like, frank_like_df, GDP_frank_pred, frank_pred_val, frank_list_pred_lb, frank_pred_val_vec)
  gc()
  
  pred_cond_mod_tau = pred_cond_thin %>%
    group_by(obs) %>%
    summarise(gauss_tau_mean = mean(gauss_tau), gauss_tau_q975 = quantile(gauss_tau, .975), gauss_tau_q025 = quantile(gauss_tau, .025),
              frank_tau_mean = mean(frank_tau), frank_tau_q975 = quantile(frank_tau, .975), frank_tau_q025 = quantile(frank_tau, .025),
              t_tau_mean = mean(t_tau), t_tau_q975 = quantile(t_tau, .975), t_tau_q025 = quantile(t_tau, .025),
              clayton_tau_mean = mean(clayton_tau), clayton_tau_q975 = quantile(clayton_tau, .975), clayton_tau_q025 = quantile(clayton_tau, .025),
              gumbel_tau_mean = mean(gumbel_tau), gumbel_tau_q975 = quantile(gumbel_tau, .975), gumbel_tau_q025 = quantile(gumbel_tau, .025))
  
  pl_tau_est <- ggplot(pred_cond_mod_tau) +
    geom_line(aes(obs, gauss_tau_mean, col = "gauss")) +
    geom_line(aes(obs, gauss_tau_q975, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, gauss_tau_q025, col = "gauss"),linetype="dotted") +
    geom_line(aes(obs, t_tau_mean, col = "t")) +
    geom_line(aes(obs, t_tau_q975, col = "t"),linetype="dotted") +
    geom_line(aes(obs, t_tau_q025, col = "t"),linetype="dotted") +
    geom_line(aes(obs, clayton_tau_mean, col = "clayton")) +
    geom_line(aes(obs, clayton_tau_q975, col = "clayton"),linetype="dotted") +
    geom_line(aes(obs, clayton_tau_q025, col = "clayton"),linetype="dotted") +
    geom_line(aes(obs, gumbel_tau_mean, col = "gumbel")) +
    geom_line(aes(obs, gumbel_tau_q975, col = "gumbel"),linetype="dotted") +
    geom_line(aes(obs, gumbel_tau_q025, col = "gumbel"),linetype="dotted") +
    geom_line(aes(obs, frank_tau_mean, col = "frank")) +
    geom_line(aes(obs, frank_tau_q975, col = "frank"),linetype="dotted") +
    geom_line(aes(obs, frank_tau_q025, col = "frank"),linetype="dotted") +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pl_tau_est
  
  save(pred_cond, pred_cond_mod_tau, pl_tau_est, file = "pred_tau_LT_adapt.Rdata")
  
  cor(U1_LT,U2_LT, method = "kendall")
  
  par(mar = c(5,5,2,1), mfrow = c(2,3))
  
  plot(U1_LT,U2_LT, xlab = "Female Literacy", ylab = "Male Literacy", main = "Observed data")
  plot(gauss_pred_U1_LT,gauss_pred_U2_LT, xlab = "Female Literacy", ylab = "Male Literacy", main = "Predicted (Gaussian)")
  plot(t_pred_U1_LT,t_pred_U2_LT, xlab = "Female Literacy", ylab = "Male Literacy", main = "Predicted (Student-t)")
  plot(clayton_pred_U1_LT,clayton_pred_U2_LT, xlab = "Female Literacy", ylab = "Male Literacy", main = "Predicted (Clayton)")
  plot(gumbel_pred_U1_LT,gumbel_pred_U2_LT, xlab = "Female Literacy", ylab = "Male Literacy", main = "Predicted (Gumbel)")
  plot(frank_pred_U1_LT,frank_pred_U2_LT, xlab = "Female Literacy", ylab = "Male Literacy", main = "Predicted (Frank)")
  # 10 7
  
  # helper for test statistics
  
  library(cramer)
  
  cram_gauss <- cramer.test(cbind(U1_LT,U2_LT), cbind(gauss_pred_U1_LT, gauss_pred_U2_LT), replicates = 10000)
  cram_t <- cramer.test(cbind(U1_LT,U2_LT), cbind(t_pred_U1_LT, t_pred_U2_LT), replicates = 10000, sim = "permutation")
  cram_clayton <- cramer.test(cbind(U1_LT,U2_LT), cbind(clayton_pred_U1_LT, clayton_pred_U2_LT), replicates = 10000)
  cram_gumbel <- cramer.test(cbind(U1_LT,U2_LT), cbind(gumbel_pred_U1_LT, gumbel_pred_U2_LT), replicates = 10000)
  cram_frank <- cramer.test(cbind(U1_LT,U2_LT), cbind(frank_pred_U1_LT, frank_pred_U2_LT), replicates = 10000)
  
  library(fasano.franceschini.test)
  
  ff_gauss <- fasano.franceschini.test(cbind(U1_LT,U2_LT), cbind(gauss_pred_U1_LT, gauss_pred_U2_LT), nPermute = 10000)
  ff_t <- fasano.franceschini.test(cbind(U1_LT,U2_LT), cbind(t_pred_U1_LT, t_pred_U2_LT), nPermute = 10000)
  ff_clayton <- fasano.franceschini.test(cbind(U1_LT,U2_LT), cbind(clayton_pred_U1_LT, clayton_pred_U2_LT), nPermute = 10000)
  ff_gumbel <- fasano.franceschini.test(cbind(U1_LT,U2_LT), cbind(gumbel_pred_U1_LT, gumbel_pred_U2_LT), nPermute = 10000)
  ff_frank <- fasano.franceschini.test(cbind(U1_LT,U2_LT), cbind(frank_pred_U1_LT, frank_pred_U2_LT), nPermute = 10000)
  
  library(ks)
  
  kde_gauss <- kde.test(cbind(U1_LT,U2_LT), cbind(gauss_pred_U1_LT, gauss_pred_U2_LT))
  kde_t <- kde.test(cbind(U1_LT,U2_LT), cbind(t_pred_U1_LT, t_pred_U2_LT))
  kde_clayton <- kde.test(cbind(U1_LT,U2_LT), cbind(clayton_pred_U1_LT, clayton_pred_U2_LT))
  kde_gumbel <- kde.test(cbind(U1_LT,U2_LT), cbind(gumbel_pred_U1_LT, gumbel_pred_U2_LT))
  kde_frank <- kde.test(cbind(U1_LT,U2_LT), cbind(frank_pred_U1_LT, frank_pred_U2_LT))
  
  p_val_summ <- rbind(c(cram_gauss$p.value, ff_gauss$p.value, kde_gauss$pvalue),
                      c(cram_t$p.value, ff_t$p.value, kde_t$pvalue),
                      c(cram_clayton$p.value, ff_clayton$p.value, kde_clayton$pvalue),
                      c(cram_gumbel$p.value, ff_gumbel$p.value, kde_gumbel$pvalue),
                      c(cram_frank$p.value, ff_frank$p.value, kde_frank$pvalue))
  
  xtable(p_val_summ)
  
  hist_true <- hist2d(U1_LT, U2_LT, nbins = c(10,10), show = FALSE)
  
  
  hist3D(
    x = hist_true$x,
    y = hist_true$y,
    z = hist_true$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.02, bty = "g", phi = 30,    # shading gives 3D effect
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
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.02, bty = "g", phi = 30,    # shading gives 3D effect
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
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.02, bty = "g", phi = 30,    # shading gives 3D effect
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
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.02, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "Clayton"
  )
  
  hist3D(
    x = hist_gumbel$x,
    y = hist_gumbel$y,
    z = hist_gumbel$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.02, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "gumbel"
  )
  
  hist3D(
    x = hist_frank$x,
    y = hist_frank$y,
    z = hist_frank$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "grey",
    theta = -45, scale = FALSE, expand = 0.02, bty = "g", phi = 30,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "Frank"
  )
  
  ################################################################################
  
}

################################################################################
# comparative plotting
################################################################################

if(F){
  wo_adapt <- new.env()
  load("LE_vs_GDP_post_process/gauss_gdp_LE_post.Rdata", envir = wo_adapt)
  gauss_woa <- as.list(wo_adapt)
  wo_adapt <- new.env()
  load("LE_vs_GDP_post_process/t_gdp_LE_post.Rdata", envir = wo_adapt)
  t_woa <- as.list(wo_adapt)
  wo_adapt <- new.env()
  load("LE_vs_GDP_post_process/clayton_gdp_LE_post.Rdata", envir = wo_adapt)
  clayton_woa <- as.list(wo_adapt)
  wo_adapt <- new.env()
  load("LE_vs_GDP_post_process/gumbel_gdp_LE_post.Rdata", envir = wo_adapt)
  gumbel_woa <- as.list(wo_adapt)
  wo_adapt <- new.env()
  load("LE_vs_GDP_post_process/frank_gdp_LE_post.Rdata", envir = wo_adapt)
  frank_woa <- as.list(wo_adapt)
  wo_adapt <- new.env()
  load("LE_vs_GDP_post_process/pred_tau_LE.Rdata", envir = wo_adapt)
  pred_woa <- as.list(wo_adapt)
  rm(wo_adapt)
  w_adapt <- new.env()
  load("LE_vs_GDP_post_process_adapt/gauss_gdp_LE_post_adapt.Rdata", envir = w_adapt)
  gauss_wa <- as.list(w_adapt)
  w_adapt <- new.env()
  load("LE_vs_GDP_post_process_adapt/t_gdp_LE_post_adapt.Rdata", envir = w_adapt)
  t_wa <- as.list(w_adapt)
  w_adapt <- new.env()
  load("LE_vs_GDP_post_process_adapt/clayton_gdp_LE_post_adapt.Rdata", envir = w_adapt)
  clayton_wa <- as.list(w_adapt)
  w_adapt <- new.env()
  load("LE_vs_GDP_post_process_adapt/gumbel_gdp_LE_post_adapt.Rdata", envir = w_adapt)
  gumbel_wa <- as.list(w_adapt)
  w_adapt <- new.env()
  load("LE_vs_GDP_post_process_adapt/frank_gdp_LE_post_adapt.Rdata", envir = w_adapt)
  frank_wa <- as.list(w_adapt)
  w_adapt <- new.env()
  load("LE_vs_GDP_post_process_adapt/pred_tau_LE_adapt.Rdata", envir = w_adapt)
  pred_wa <- as.list(w_adapt)
  rm(w_adapt)
  
  
  
  (gauss_woa$gauss_pl_like + labs(title="Gaussian (without adaption)") + gauss_wa$gauss_pl_like + labs(title="Gaussian (with adaption)")) / 
    (t_woa$t_pl_like + labs(title="Student-t (without adaption)") + t_wa$t_pl_like + labs(title="Student-t (with adaption)")) /
    (clayton_woa$clayton_pl_like + labs(title="Clayton (without adaption)") + clayton_wa$clayton_pl_like + labs(title="Clayton (with adaption)")) /
    (gumbel_woa$gumbel_pl_like + labs(title="Gumbel (without adaption)") + gumbel_wa$gumbel_pl_like + labs(title="Gumbel (with adaption)")) / 
    (frank_woa$frank_pl_like + labs(title="Frank (without adaption)") + frank_wa$frank_pl_like + labs(title="Frank (with adaption)"))
  
  pred_woa$pl_tau_est + labs(title="Without adaption") + xlab("Scaled GDP") + ylim(-0.6,1) + pred_wa$pl_tau_est + labs(title="With adaption") + ylim(-0.6,1) + xlab("Scaled GDP")
  
}