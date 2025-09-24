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
  
  load("gauss_gdp_LE.Rdata")
  
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
  
  GDP_gauss_pred <- BiCopSim(N = nrow(GDP), family = 1, par = link_gauss(colMeans(gauss_pred_val)))
  
  gauss_pred_U1_LE = GDP_gauss_pred[,1]
  gauss_pred_U2_LE = GDP_gauss_pred[,2]
  
  hist_gauss <- hist2d(gauss_pred_U1_LE, gauss_pred_U2_LE, nbins = c(10,10), show = FALSE)
  
  save(gauss_pred_U1_LE, gauss_pred_U2_LE, hist_gauss, gauss_like_val, gauss_pl_like, gauss_like_df, gauss_pred_val_vec, file = "gauss_gdp_LE_post.Rdata")
  rm(gauss_GDP_LE, gauss_pred_U1_LE, gauss_pred_U2_LE, hist_gauss, gauss_like_val, gauss_pl_like, gauss_like_df, GDP_gauss_pred, gauss_pred_val, gauss_list_pred_lb, gauss_pred_val_vec)
  gc()
  
  load("frank_gdp_LE.Rdata")
  
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
  
  GDP_frank_pred <- BiCopSim(N = nrow(GDP), family = 1, par = link_frank(colMeans(frank_pred_val)))
  
  frank_pred_U1_LE = GDP_frank_pred[,1]
  frank_pred_U2_LE = GDP_frank_pred[,2]
  
  hist_frank <- hist2d(frank_pred_U1_LE, frank_pred_U2_LE, nbins = c(10,10), show = FALSE)
  
  save(frank_pred_U1_LE, frank_pred_U2_LE, hist_frank, frank_like_val, frank_pl_like, frank_like_df, frank_pred_val_vec, file = "frank_gdp_LE_post.Rdata")
  rm(frank_GDP_LE, frank_pred_U1_LE, frank_pred_U2_LE, hist_frank, frank_like_val, frank_pl_like, frank_like_df, GDP_frank_pred, frank_pred_val, frank_list_pred_lb, frank_pred_val_vec)
  gc()
  
  load("t_gdp_LE.Rdata")
  
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
  
  GDP_t_pred <- BiCopSim(N = nrow(GDP), family = 1, par = link_t(colMeans(t_pred_val)))
  
  t_pred_U1_LE = GDP_t_pred[,1]
  t_pred_U2_LE = GDP_t_pred[,2]
  
  hist_t <- hist2d(t_pred_U1_LE, t_pred_U2_LE, nbins = c(10,10), show = FALSE)
  
  save(t_pred_U1_LE, t_pred_U2_LE, hist_t, t_like_val, t_pl_like, t_like_df, t_pred_val_vec, file = "t_gdp_LE_post.Rdata")
  rm(t_GDP_LE, t_pred_U1_LE, t_pred_U2_LE, hist_t, t_like_val, t_pl_like, t_like_df, GDP_t_pred, t_pred_val, t_list_pred_lb, t_pred_val_vec)
  gc()
  
  
  
  load("clayton_gdp_LE.Rdata")
  
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
  
  GDP_clayton_pred <- BiCopSim(N = nrow(GDP), family = 1, par = link_clayton(colMeans(clayton_pred_val)))
  
  clayton_pred_U1_LE = GDP_clayton_pred[,1]
  clayton_pred_U2_LE = GDP_clayton_pred[,2]
  
  hist_clayton <- hist2d(clayton_pred_U1_LE, clayton_pred_U2_LE, nbins = c(10,10), show = FALSE)
  
  save(clayton_pred_U1_LE, clayton_pred_U2_LE, hist_clayton, clayton_like_val, clayton_pl_like, clayton_like_df, clayton_pred_val_vec, file = "clayton_gdp_LE_post.Rdata")
  rm(clayton_GDP_LE, clayton_pred_U1_LE, clayton_pred_U2_LE, hist_clayton, clayton_like_val, clayton_pl_like, clayton_like_df, GDP_clayton_pred, clayton_pred_val, clayton_list_pred_lb, clayton_pred_val_vec)
  gc()
  
  load("gumbel_gdp_LE.Rdata")
  
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
  
  GDP_gumbel_pred <- BiCopSim(N = nrow(GDP), family = 1, par = link_gumbel(colMeans(gumbel_pred_val)))
  
  gumbel_pred_U1_LE = GDP_gumbel_pred[,1]
  gumbel_pred_U2_LE = GDP_gumbel_pred[,2]
  
  hist_gumbel <- hist2d(gumbel_pred_U1_LE, gumbel_pred_U2_LE, nbins = c(10,10), show = FALSE)
  
  save(gumbel_pred_U1_LE, gumbel_pred_U2_LE, hist_gumbel, gumbel_like_val, gumbel_pl_like, gumbel_like_df, gumbel_pred_val_vec, file = "gumbel_gdp_LE_post.Rdata")
  rm(gumbel_GDP_LE, gumbel_pred_U1_LE, gumbel_pred_U2_LE, hist_gumbel, gumbel_like_val, gumbel_pl_like, gumbel_like_df, GDP_gumbel_pred, gumbel_pred_val, gumbel_list_pred_lb, gumbel_pred_val_vec)
  gc()
  
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
  
  save(pred_cond, pred_cond_mod_tau, pl_tau_est, file = "pred_tau_LE.Rdata")
  
  cor(U1_LE,U2_LE, method = "kendall")
  cor(U1_LE,U2_LE, method = "kendall")
  
  plot(U1_LE,U2_LE)
  plot(gauss_pred_U1_LE,gauss_pred_U2_LE)
  
  plot(U1_LE,U2_LE)
  plot(t_pred_U1_LE,t_pred_U2_LE)
  
  plot(U1_LE,U2_LE)
  plot(clayton_pred_U1_LE,clayton_pred_U2_LE)

  plot(U1_LE,U2_LE)
  plot(gumbel_pred_U1_LE,gumbel_pred_U2_LE)
  
  plot(U1_LE,U2_LE)
  plot(frank_pred_U1_LE,frank_pred_U2_LE)
  
  hist_true <- hist2d(U1_LE, U2_LE, nbins = c(10,10), show = FALSE)
  
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
