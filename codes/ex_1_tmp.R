library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

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
n.iter_par <- 50000
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

# load("50k_mcmc/gauss_gdp_LE_tree_10.Rdata")
# 
# gauss_pred_val <- do.call(rbind,lapply(1:length(gauss_GDP_LE$trees), \(idx) BART_calculate_pred(gauss_GDP_LE$trees[[idx]], GDP)))
# 
# gauss_like_df <-data.frame("nn" = apply(gauss_pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1_LE, U2_LE)))
# gauss_like_df$idx <- rep(1:n.iter_par, n.chain_par)
# gauss_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
# 
# gauss_pl_like <- gauss_like_df %>%
#   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
#   geom_line() +
#   labs(
#     x = "Iteration",
#     y = "Log-likelihood"
#   ) +
#   guides(color = "none") +
#   theme_minimal()
# 
# save(gauss_pl_like, file = "plot_gauss_LE.Rdata")
# 
# rm(gauss_pred_val, gauss_like_df, gauss_pl_like, gauss_GDP_LE)
# 
# gc()
# rstudioapi::restartSession()
# 
# load("50k_mcmc/gauss_gdp_LE_tree_10_adapt.Rdata")
# 
# gauss_pred_val <- do.call(rbind,lapply(1:length(gauss_GDP_LE$trees), \(idx) BART_calculate_pred(gauss_GDP_LE$trees[[idx]], GDP)))
# 
# gauss_like_df <-data.frame("nn" = apply(gauss_pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1_LE, U2_LE)))
# gauss_like_df$idx <- rep(1:n.iter_par, n.chain_par)
# gauss_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
# 
# gauss_pl_like <- gauss_like_df %>%
#   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
#   geom_line() +
#   labs(
#     x = "Iteration",
#     y = "Log-likelihood"
#   ) +
#   guides(color = "none") +
#   theme_minimal()
# 
# save(gauss_pl_like, file = "plot_gauss_LE_adapt.Rdata")
# 
# rm(gauss_pred_val, gauss_like_df, gauss_pl_like, gauss_GDP_LE)
# 
# gc()
# rstudioapi::restartSession()
# 
# load("50k_mcmc/gauss_gdp_LT_tree_10.Rdata")
# 
# gauss_pred_val <- do.call(rbind,lapply(1:length(gauss_GDP_LT$trees), \(idx) BART_calculate_pred(gauss_GDP_LT$trees[[idx]], GDP)))
# 
# gauss_like_df <-data.frame("nn" = apply(gauss_pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1_LT, U2_LT)))
# gauss_like_df$idx <- rep(1:n.iter_par, n.chain_par)
# gauss_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
# 
# gauss_pl_like <- gauss_like_df %>%
#   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
#   geom_line() +
#   labs(
#     x = "Iteration",
#     y = "Log-likelihood"
#   ) +
#   guides(color = "none") +
#   theme_minimal()
# 
# save(gauss_pl_like, file = "plot_gauss_LT.Rdata")
# 
# rm(gauss_pred_val, gauss_like_df, gauss_pl_like, gauss_GDP_LT)
# 
# gc()
# rstudioapi::restartSession()
# 
# load("50k_mcmc/gauss_gdp_LT_tree_10_adapt.Rdata")
# 
# gauss_pred_val <- do.call(rbind,lapply(1:length(gauss_GDP_LT$trees), \(idx) BART_calculate_pred(gauss_GDP_LT$trees[[idx]], GDP)))
# 
# gauss_like_df <-data.frame("nn" = apply(gauss_pred_val, 1, function(x)loglik_gauss(link_gauss(x), U1_LT, U2_LT)))
# gauss_like_df$idx <- rep(1:n.iter_par, n.chain_par)
# gauss_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
# 
# gauss_pl_like <- gauss_like_df %>%
#   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
#   geom_line() +
#   labs(
#     x = "Iteration",
#     y = "Log-likelihood"
#   ) +
#   guides(color = "none") +
#   theme_minimal()
# 
# save(gauss_pl_like, file = "plot_gauss_LT_adapt.Rdata")
# 
# rm(gauss_pred_val, gauss_like_df, gauss_pl_like, gauss_GDP_LT)
# 
# gc()
# rstudioapi::restartSession()
# 
# load("50k_mcmc/t_gdp_LE_tree_10.Rdata")
# 
# t_pred_val <- do.call(rbind,lapply(1:length(t_GDP_LE$trees), \(idx) BART_calculate_pred(t_GDP_LE$trees[[idx]], GDP)))
# 
# t_like_df <-data.frame("nn" = apply(t_pred_val, 1, function(x)loglik_t(link_t(x), U1_LE, U2_LE)))
# t_like_df$idx <- rep(1:n.iter_par, n.chain_par)
# t_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
# 
# t_pl_like <- t_like_df %>%
#   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
#   geom_line() +
#   labs(
#     x = "Iteration",
#     y = "Log-likelihood"
#   ) +
#   guides(color = "none") +
#   theme_minimal()
# 
# save(t_pl_like, file = "plot_t_LE.Rdata")
# 
# rm(t_pred_val, t_like_df, t_pl_like, t_GDP_LE)
# 
# gc()
# rstudioapi::restartSession()
# 
# load("50k_mcmc/t_gdp_LE_tree_10_adapt.Rdata")
# 
# t_pred_val <- do.call(rbind,lapply(1:length(t_GDP_LE$trees), \(idx) BART_calculate_pred(t_GDP_LE$trees[[idx]], GDP)))
# 
# t_like_df <-data.frame("nn" = apply(t_pred_val, 1, function(x)loglik_t(link_t(x), U1_LE, U2_LE)))
# t_like_df$idx <- rep(1:n.iter_par, n.chain_par)
# t_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
# 
# t_pl_like <- t_like_df %>%
#   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
#   geom_line() +
#   labs(
#     x = "Iteration",
#     y = "Log-likelihood"
#   ) +
#   guides(color = "none") +
#   theme_minimal()
# 
# save(t_pl_like, file = "plot_t_LE_adapt.Rdata")
# 
# rm(t_pred_val, t_like_df, t_pl_like, t_GDP_LE)
# 
# gc()
# rstudioapi::restartSession()
# 
# load("50k_mcmc/t_gdp_LT_tree_10.Rdata")
# 
# t_pred_val <- do.call(rbind,lapply(1:length(t_GDP_LT$trees), \(idx) BART_calculate_pred(t_GDP_LT$trees[[idx]], GDP)))
# 
# t_like_df <-data.frame("nn" = apply(t_pred_val, 1, function(x)loglik_t(link_t(x), U1_LT, U2_LT)))
# t_like_df$idx <- rep(1:n.iter_par, n.chain_par)
# t_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
# 
# t_pl_like <- t_like_df %>%
#   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
#   geom_line() +
#   labs(
#     x = "Iteration",
#     y = "Log-likelihood"
#   ) +
#   guides(color = "none") +
#   theme_minimal()
# 
# save(t_pl_like, file = "plot_t_LT.Rdata")
# 
# rm(t_pred_val, t_like_df, t_pl_like, t_GDP_LT)
# 
# gc()
# rstudioapi::restartSession()
# 
# load("50k_mcmc/t_gdp_LT_tree_10_adapt.Rdata")
# 
# t_pred_val <- do.call(rbind,lapply(1:length(t_GDP_LT$trees), \(idx) BART_calculate_pred(t_GDP_LT$trees[[idx]], GDP)))
# 
# t_like_df <-data.frame("nn" = apply(t_pred_val, 1, function(x)loglik_t(link_t(x), U1_LT, U2_LT)))
# t_like_df$idx <- rep(1:n.iter_par, n.chain_par)
# t_like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
# 
# t_pl_like <- t_like_df %>%
#   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
#   geom_line() +
#   labs(
#     x = "Iteration",
#     y = "Log-likelihood"
#   ) +
#   guides(color = "none") +
#   theme_minimal()
# 
# save(t_pl_like, file = "plot_t_LT_adapt.Rdata")
# 
# rm(t_pred_val, t_like_df, t_pl_like, t_GDP_LT)
# 
# gc()
# rstudioapi::restartSession()

gauss_LT <- get(load("plot_gauss_LT.Rdata"))
gauss_LT_adapt <- get(load("plot_gauss_LT_adapt.Rdata"))
t_LT <- get(load("plot_t_LT.Rdata"))
t_LT_adapt <- get(load("plot_t_LT_adapt.Rdata"))

gauss_LE <- get(load("plot_gauss_LE.Rdata"))
gauss_LE_adapt <- get(load("plot_gauss_LE_adapt.Rdata"))
t_LE <- get(load("plot_t_LE.Rdata"))
t_LE_adapt <- get(load("plot_t_LE_adapt.Rdata"))


(gauss_LE + labs(title="LE-Gaussian (without adaption)") + ylim(0,400) + gauss_LE_adapt + labs(title="LE-Gaussian (with adaption)") + ylim(0,400)) / 
  (t_LE + labs(title="LE-Student-t (without adaption)") + ylim(50,400) + t_LE_adapt + labs(title="LE-Student-t (with adaption)") + ylim(50,400)) /
  (gauss_LT + labs(title="LT-Gaussian (without adaption)") + ylim(50,400) + gauss_LT_adapt + labs(title="LT-Gaussian (with adaption)") + ylim(50,400)) /
  (t_LT + labs(title="LT-Student-t (without adaption)") + ylim(150,450) + t_LT_adapt + labs(title="LT-Student-t (with adaption)") + ylim(150,450)) 
