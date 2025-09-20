gauss <- get(load("gauss_mcmc_1_tree_1_plot.Rdata"))
gauss_stat <- get(load("gauss_mcmc_1_tree_1_stat.Rdata"))

gauss_adapt <- get(load("gauss_mcmc_1_tree_1_plot_adapt.Rdata"))
gauss_stat_adapt <- get(load("gauss_mcmc_1_tree_1_stat_adapt.Rdata"))

gauss_like_true <- loglik_gauss((sin(tau_true_1 * pi/2)), copula_uu_gauss_1[,1] , copula_uu_gauss_1[,2])
# gauss_like_true <- loglik_gauss((sin(tau_true_2 * pi/2)), copula_uu_gauss_2[,1] , copula_uu_gauss_2[,2])


t <- get(load("t_mcmc_1_tree_1_plot.Rdata"))
t_stat <- get(load("t_mcmc_1_tree_1_stat.Rdata"))

t_adapt <- get(load("t_mcmc_1_tree_1_plot_adapt.Rdata"))
t_stat_adapt <- get(load("t_mcmc_1_tree_1_stat_adapt.Rdata"))

t_like_true <- loglik_t((sin(tau_true_1 * pi/2)), copula_uu_t_1[,1] , copula_uu_t_1[,2])
# t_like_true <- loglik_t((sin(tau_true_2 * pi/2)), copula_uu_t_2[,1] , copula_uu_t_2[,2])


clayton <- get(load("clayton_mcmc_1_tree_1_plot.Rdata"))
clayton_stat <- get(load("clayton_mcmc_1_tree_1_stat.Rdata"))

clayton_adapt <- get(load("clayton_mcmc_1_tree_1_plot_adapt.Rdata"))
clayton_stat_adapt <- get(load("clayton_mcmc_1_tree_1_stat_adapt.Rdata"))

clayton_like_true <- loglik_clayton((2*tau_true_1 / (1-tau_true_1)), copula_uu_clayton_1[,1] , copula_uu_clayton_1[,2])
# clayton_like_true <- loglik_clayton((2*tau_true_2 / (1-tau_true_2)), copula_uu_clayton_2[,1] , copula_uu_clayton_2[,2])


gumbel <- get(load("gumbel_mcmc_1_tree_1_plot.Rdata"))
gumbel_stat <- get(load("gumbel_mcmc_1_tree_1_stat.Rdata"))

gumbel_adapt <- get(load("gumbel_mcmc_1_tree_1_plot_adapt.Rdata"))
gumbel_stat_adapt <- get(load("gumbel_mcmc_1_tree_1_stat_adapt.Rdata"))

gumbel_like_true <- loglik_gumbel((1 / (1-tau_true_1)), copula_uu_gumbel_1[,1] , copula_uu_gumbel_1[,2])
# gumbel_like_true <- loglik_gumbel((1 / (1-tau_true_2)), copula_uu_gumbel_2[,1] , copula_uu_gumbel_2[,2])


frank <- get(load("frank_mcmc_1_tree_1_plot.Rdata"))
frank_stat <- get(load("frank_mcmc_1_tree_1_stat.Rdata"))

frank_adapt <- get(load("frank_mcmc_1_tree_1_plot_adapt.Rdata"))
frank_stat_adapt <- get(load("frank_mcmc_1_tree_1_stat_adapt.Rdata"))


frank_like_true <- loglik_frank(BiCopTau2Par(5,tau_true_1), copula_uu_frank_1[,1] , copula_uu_frank_1[,2])
# frank_like_true <- loglik_frank(BiCopTau2Par(5,tau_true_2), copula_uu_frank_2[,1] , copula_uu_frank_2[,2])

c(gauss_like_true, t_like_true, clayton_like_true, gumbel_like_true, frank_like_true)

summ_stat_acc <- rbind(c(as.vector(gauss_stat$pred),as.vector(gauss_stat_adapt$pred)),
                       c(as.vector(t_stat$pred),as.vector(t_stat_adapt$pred)),
                       c(as.vector(clayton_stat$pred),as.vector(clayton_stat_adapt$pred)),
                       c(as.vector(gumbel_stat$pred),as.vector(gumbel_stat_adapt$pred)),
                       c(as.vector(frank_stat$pred),as.vector(frank_stat_adapt$pred)))

xtable(summ_stat_acc, digits = 4)

summ_stat_tree <- rbind(c(as.vector(gauss_stat$tree),as.vector(gauss_stat_adapt$tree)),
                        c(as.vector(t_stat$tree),as.vector(t_stat_adapt$tree)),
                        c(as.vector(clayton_stat$tree),as.vector(clayton_stat_adapt$tree)),
                        c(as.vector(gumbel_stat$tree),as.vector(gumbel_stat_adapt$tree)),
                        c(as.vector(frank_stat$tree),as.vector(frank_stat_adapt$tree)))

xtable(summ_stat_tree, digits = 2)

(gauss$like + labs(title="Gaussian (without adaption)") + geom_hline(yintercept = gauss_like_true, linetype = 2) + gauss_adapt$like + labs(title="Gaussian (with adaption)") + geom_hline(yintercept = gauss_like_true, linetype = 2)) / 
  (t$like + labs(title="Student-t (without adaption)") + geom_hline(yintercept = t_like_true, linetype = 2) + t_adapt$like + labs(title="Student-t (with adaption)") + geom_hline(yintercept = t_like_true, linetype = 2)) /
  (clayton$like + labs(title="Clayton (without adaption)")+ geom_hline(yintercept = clayton_like_true, linetype = 2) + clayton_adapt$like + labs(title="Clayton (with adaption)")+ geom_hline(yintercept = clayton_like_true, linetype = 2)) /
  (gumbel$like + labs(title="Gumbel (without adaption)")+ geom_hline(yintercept = gumbel_like_true, linetype = 2) + gumbel_adapt$like + labs(title="Gumbel (with adaption)")+ geom_hline(yintercept = gumbel_like_true, linetype = 2)) /
  (frank$like + labs(title="Frank (without adaption)")+ geom_hline(yintercept = frank_like_true, linetype = 2) + frank_adapt$like + labs(title="Frank (with adaption)")+ geom_hline(yintercept = frank_like_true, linetype = 2))

(gauss$pred + labs(title="Gaussian (without adaption)") + ylim(-0.2,1) + gauss_adapt$pred + labs(title="Gaussian (with adaption)") + ylim(-0.2,1)) / 
  (t$pred + labs(title="Student-t (without adaption)") + ylim(-0.2,1) + t_adapt$pred + labs(title="Student-t (with adaption)") + ylim(-0.2,1)) /
  (clayton$pred + labs(title="Clayton (without adaption)") + ylim(0,1) + clayton_adapt$pred + labs(title="Clayton (with adaption)") + ylim(0,1)) /
  (gumbel$pred + labs(title="Gumbel (without adaption)") + ylim(0,1) + gumbel_adapt$pred + labs(title="Gumbel (with adaption)") + ylim(0,1)) /
  (frank$pred + labs(title="Frank (without adaption)") + ylim(-0.1,1) + frank_adapt$pred + labs(title="Frank (with adaption)") + ylim(-0.1,1))

# (gauss$nterm + labs(title="Gaussian (without adaption)") + geom_hline(yintercept = 3, linetype = 2) + ylim(0, 8) + gauss_adapt$nterm + labs(title="Gaussian (with adaption)") + geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) / 
#   (t$nterm + labs(title="Student-t (without adaption)") + geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + t_adapt$nterm + labs(title="Student-t (with adaption)") + geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) /
#   (clayton$nterm + labs(title="Clayton (without adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + clayton_adapt$nterm + labs(title="Clayton (with adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) /
#   (gumbel$nterm + labs(title="Gumbel (without adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + gumbel_adapt$nterm + labs(title="Gumbel (with adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) /
#   (frank$nterm + labs(title="Frank (without adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + frank_adapt$nterm + labs(title="Frank (with adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8))
# 
# (gauss$depth + labs(title="Gaussian (without adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + gauss_adapt$depth + labs(title="Gaussian (with adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) / 
#   (t$depth + labs(title="Student-t (without adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + t_adapt$depth + labs(title="Student-t (with adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) /
#   (clayton$depth + labs(title="Clayton (without adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + clayton_adapt$depth + labs(title="Clayton (with adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) /
#   (gumbel$depth + labs(title="Gumbel (without adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + gumbel_adapt$depth + labs(title="Gumbel (with adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) /
#   (frank$depth + labs(title="Frank (without adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + frank_adapt$depth + labs(title="Frank (with adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6))
