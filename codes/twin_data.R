# install.packages("misty")

library(misty)
library(dplyr)

var_twin1 <- c("V10005", "V10006", "V10007", "V10008", "V10009")
var_twin2 <- c("V20005", "V20006", "V20007", "V20008", "V20009")
var_pars <- c("V30009", "V30010", "V30011")

all_data <- read.sav("dataverse_files/nmgen.sav")

twin1_data <- all_data[,colnames(all_data)%in%var_twin1]
twin2_data <- all_data[,colnames(all_data)%in%var_twin2]
pars_data <- all_data[,colnames(all_data) %in%var_pars]

red_all_data <- cbind(rowSums(twin1_data), rowSums(twin2_data), pars_data)

red_all_data <- na.omit(red_all_data)

colnames(red_all_data) <- c("T1_score", "T2_score", "M_ed", "F_ed", "Income")

U1 = ecdf(red_all_data$T1_score)(red_all_data$T1_score)
U2 = ecdf(red_all_data$T2_score)(red_all_data$T2_score)

M_ed <- as.data.frame((red_all_data$M_ed - min(red_all_data$M_ed))/(max(red_all_data$M_ed) - min(red_all_data$M_ed)))
M_ed <- as.matrix(M_ed)
rownames(M_ed) <- 1:nrow(M_ed)

plot(U1,U2)

n.chain_par <- 1
n.iter_par <- 6000
n.born.out.par <- 1000
n.thin <- 1
incl.split_par <- TRUE
cont.unif_par <- TRUE
moves.prob_par <- c(0.4, 0.4, 0.1, 0.1)
##########################################################

gauss_M_ed_tree_1 <- MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                 n.tree = 1,
            X = M_ed,
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


model_twin <- gauss_M_ed_tree_1

list_pred_lb <- lapply(1:length(model_twin$trees), \(idx) BART_calculate_pred(model_twin$trees[[idx]], M_ed))

pred_val = do.call(rbind,list_pred_lb)

n.thin <- 1
n.iter_par <- 6000
n.born.out.par <- 1000

pred_val_vec = as.vector(pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])

pred_obs = rep(M_ed, each = (n.chain_par * (n.iter_par - n.born.out.par)))

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

M_ed_gauss_pred <- BiCopSim(N = nrow(M_ed), family = 1, par = link_gauss(colMeans(pred_val)))

pred_U1 = M_ed_gauss_pred[,1]
pred_U2 = M_ed_gauss_pred[,2]

plot(U1,U2)
plot(pred_U1,pred_U2)

