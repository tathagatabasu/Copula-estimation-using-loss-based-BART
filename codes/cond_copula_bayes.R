# codes and packages

#install.packages(c("copula", "MASS", "coda"))
source('code/import_functions.R')
library(data.tree)
library(ggplot2)
library(ggpubr)
library(copula)
library(MASS)   # For multivariate normal functions
library(coda)   # For MCMC diagnostics

# data generation

set.seed(123)

# generate predictor
n <- 500
x_pred <- rnorm(n)

# Define true copula parameter
calib_true <- 1+2*x_pred^2

rho_true <- 2 /(abs(calib_true)+1)-1

# Define the Gaussian copula
cop <- sapply(rho_true, function(x)normalCopula(param = x, dim = 2))

# Generate uniform data from the copula
u <- lapply(cop, function(x)rCopula(1, x))
u <- matrix(unlist(u), ncol = 2, byrow = T)

# Convert to normal marginals
y1 <- qnorm(u[,1])
y2 <- qnorm(u[,2])

# empirical cdf function
ecdf_y1 = ecdf(y1)
ecdf_y2 = ecdf(y2)

pseudo_u1 = sapply(y1, ecdf_y1)
pseudo_u2 = sapply(y2, ecdf_y2)

triweight_ker <- function(u){
  ifelse(abs(u)<=1, 35/32*(1-u^2)^3, 0)
}

NW_ker_reg <- function(x, obs_x, band = 0.5){
  factors_NW <- sapply(obs_x, function(y)triweight_ker((y-x)/band))
  
  return(factors_NW/sum(factors_NW))
}

rho.cond=function(uu,x.vec, xgrid)
{
  uu1=uu[,1]
  uu2=uu[,2]
  estim=c()
  
  for(count in 1:length(xgrid))
  {
    # Weights
    ww=NW_ker_reg(x.vec,xgrid[count])
    
    # Compute uu.hat
    uu1.hat=c()
    uu2.hat=c()
    
    for(i in 1:nrow(uu))
    {
      uu1.hat[i]=sum(ww*(uu1<=uu1[i]))
      uu2.hat[i]=sum(ww*(uu2<=uu2[i]))
    }
    
    estim[count]=12*sum(ww*(1-uu1.hat)*(1-uu2.hat))-3
  }
  
  # Return
  return(estim)
}

sample_rho = rho.cond(cbind(y1,y2), x_pred, x_pred)

# response

theta_link_obs = 1/2*log((1+sample_rho)/(1-sample_rho))


# Store in a dataframe
data <- data.frame(d){
  
}

# define likelihood function
# import data
bcancer.df <- read.table(file = 'data/breast_cancer_wisconsin_original/breast-cancer-wisconsin.data',
                         sep = ',')
# description
#1. Sample code number            id number
#2. Clump Thickness               1 - 10
#3. Uniformity of Cell Size       1 - 10
#4. Uniformity of Cell Shape      1 - 10
#5. Marginal Adhesion             1 - 10
#6. Single Epithelial Cell Size   1 - 10
#7. Bare Nuclei                   1 - 10
#8. Bland Chromatin               1 - 10
#9. Normal Nucleoli               1 - 10
#10. Mitoses                       1 - 10
#11. Class:                        (2 for benign, 4 for malignant)

colnames(bcancer.df) <- c('ID', 
                          'clump.thick', 
                          'cell.size', 
                          'cell.shape',
                          'marg.adhesion',
                          'se.cell.size',
                          'bar.nuclei',
                          'bland.chromatin',
                          'normal.nucleoli',
                          'mitosis',
                          'cancer.class')

head(bcancer.df)
# create predictor df 
X_pred <- bcancer.df[,c('clump.thick', 
                        'cell.size', 
                        'cell.shape',
                        'marg.adhesion',
                        'se.cell.size',
                        'bar.nuclei',
                        'bland.chromatin',
                        'normal.nucleoli',
                        'mitosis')]
summary(X_pred)
# trasnform in numeric bar.nuclei
X_pred$bar.nuclei <- as.numeric(X_pred$bar.nuclei)
# remove NA 
idx.na <- is.na(X_pred$bar.nuclei)
X_pred <- X_pred[!idx.na,]
rownames(X_pred) <- 1:nrow(X_pred)
# create variable 1 if malignant and 0 if benign
Y_mal <- bcancer.df$cancer.class == 4
Y_mal <- Y_mal[!idx.na]

# create lists representing different priors - 3 omega values x 2 beta values 
lb.prior.small.hb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.small, 1.5))
lb.prior.small.lb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.small, 0.5))
lb.prior.medium.hb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.medium, 1.5))
lb.prior.medium.lb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.medium, 0.5))
lb.prior.large.hb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.large, 1.5))
lb.prior.large.lb <- list(fun = joint.prior.new.tree, param = c(omeg.chip.large, 0.5))


# normalise predictors 
X_pred.norm <- as.data.frame(apply(X_pred, 2, \(x) (x - min(x))/(max(x) - min(x))))
X_pred.norm <- as.matrix(X_pred.norm)
rownames(X_pred.norm) <- 1:nrow(X_pred.norm)
# calculate correlations

cor(X_pred.norm, Y_mal)




n.chain_par <- 100
n.iter_par <- 500
incl.split_par <- FALSE
cont.unif_par <- TRUE
moves.prob_par <- c(0.4, 0.4, 0.1, 0.1)

#############
## DEFUALT ##
#############

lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.56, 0.62))
mcmc_lb.def <- multichain_MCMC_binary(n.iter = n.iter_par,
                                      n.chain = n.chain_par,
                                      X = X_pred.norm,
                                      Y = Y_mal,
                                      alpha.prior = 1,
                                      beta.prior = 1,
                                      prior_list = lb.prior.def,
                                      include.split = incl.split_par,
                                      cont.unif = cont.unif_par,
                                      moves.prob = c(0.4, 0.4, 0.1, 0.1))

###########
## SMALL ##
###########
# - Lossbased - beta = 1.5
mcmc_lb.small.hb <- multichain_MCMC_binary(n.iter = n.iter_par,
                                           n.chain = n.chain_par,
                                           X = X_pred.norm,
                                           Y = Y_mal,
                                           alpha.prior = 1,
                                           beta.prior = 1,
                                           prior_list = lb.prior.small.hb,
                                           include.split = incl.split_par,
                                           cont.unif = cont.unif_par,
                                           moves.prob = c(0.4, 0.4, 0.1, 0.1))
# - Lossbased - beta = 0.5
mcmc_lb.small.lb <- multichain_MCMC_binary(n.iter = n.iter_par,
                                           n.chain = n.chain_par,
                                           X = X_pred.norm,
                                           Y = Y_mal,
                                           alpha.prior = 1,
                                           beta.prior = 1,
                                           prior_list = lb.prior.small.lb,
                                           include.split = incl.split_par,
                                           cont.unif = cont.unif_par,
                                           moves.prob = c(0.4, 0.4, 0.1, 0.1))

############
## MEDIUM ##
############
# - Lossbased - beta = 1.5
mcmc_lb.medium.hb <- multichain_MCMC_binary(n.iter = n.iter_par,
                                            n.chain = n.chain_par,
                                            X = X_pred.norm,
                                            Y = Y_mal,
                                            alpha.prior = 1,
                                            beta.prior = 1,
                                            prior_list = lb.prior.medium.hb,
                                            include.split = incl.split_par,
                                            cont.unif = cont.unif_par,
                                            moves.prob = c(0.4, 0.4, 0.1, 0.1))
# - Lossbased - beta = 0.5
mcmc_lb.medium.lb <- multichain_MCMC_binary(n.iter = n.iter_par,
                                            n.chain = n.chain_par,
                                            X = X_pred.norm,
                                            Y = Y_mal,
                                            alpha.prior = 1,
                                            beta.prior = 1,
                                            include.split = incl.split_par,
                                            cont.unif = cont.unif_par,
                                            prior_list = lb.prior.medium.lb,
                                            moves.prob = c(0.4, 0.4, 0.1, 0.1))
####################
## DEFAULT MODELS ##
####################

model.list.def <- list(
  mcmc_lb.def,
mcmc_lb.medium.lb,
mcmc_lb.medium.hb,
mcmc_lb.small.lb,
mcmc_lb.small.hb)

names(model.list.def) <- c(
  'LB - default',
  'LB - o = 0.30, g = 0.5',
'LB - o = 0.30, g = 1.5',
'LB - o = 0.42, g = 0.5',
'LB - o = 0.42, g = 1.5')


# extract depth, number of terminal nodes, missing rate and loglik of all the trees
depth.df <- apply_fun_models(fun_ = get_depth, 
                             mcmc.list = model.list.def,
                             born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)
nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes, 
                             mcmc.list = model.list.def, 
                             born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)




miss.rate.df <- apply_fun_models(fun_ = \(x) rapid.miss.rate(tree_top = x, X = X_pred.norm, Y = Y_mal),
                                 mcmc.list = model.list.def,
                                 born.out.pc = 250, 
                                 n.chain = n.chain_par, 
                                 sample.pc = n.iter_par)


like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_binary(tree_top = x, Y = Y_mal, 
                                                            X = X_pred.norm), 
                            mcmc.list = model.list.def, 
                            born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)




df.sum.def <- data.frame(tree = nterm.df$x,
                         panel.name = nterm.df$panel.name,
                         loglik = like.df$y,
                         nterm = nterm.df$y,
                         miss.rate = miss.rate.df$y)


df.sum.def$nterm.cut <- cut(df.sum.def$nterm, breaks = c(3,7,10,13,25))

# create data.frame for prior

df.nl.prior <- rbind(data.frame(x = 1:20, 
                                y = prior.nterm(1:20, 1.56),
                                panel.name = 'LB - default'))

# create df for post quantities
df.split <- split(df.sum.def, df.sum.def$panel.name)
post.nl.mean <- vapply(df.split, \(x) mean(x$nterm), 0)
post.nl.q0.025 <- vapply(df.split, \(x) quantile(x$nterm, 0.025), 0)
post.nl.q0.975 <- vapply(df.split, \(x) quantile(x$nterm, 0.975), 0)
df.post.nl.sum <- data.frame(mean = post.nl.mean, 
                             low = post.nl.q0.025, 
                             up = post.nl.q0.975,
                             panel.name = names(df.split))



pl.nl.post <- ggplot(df.sum.def, aes(x = nterm, y = after_stat(density))) + 
  geom_histogram(binwidth = 1, color = 'black', fill = 'white') + 
  geom_vline(xintercept = seq(1,21,by = 2), color = 'grey', alpha = 0.3) + 
  geom_line(data = df.nl.prior, aes(x,y), colour = 'deepskyblue', linewidth = 0.4) +
  geom_point(data = df.nl.prior, aes(x,y), colour = 'deepskyblue', size = 0.4) + 
  geom_vline(data = df.post.nl.sum, aes(xintercept = mean)) +
  geom_vline(data = df.post.nl.sum, aes(xintercept = low), linetype = 2) +
  geom_vline(data = df.post.nl.sum, aes(xintercept = up), linetype = 2) +
  facet_wrap(facets = ~panel.name, ncol = 1) + 
  theme_classic() + 
  xlab('Number of terminal nodes') + 
  scale_x_continuous(breaks = seq(1,21,by = 2))
pl.nl.post

pl.loglik.trace <- ggplot(df.sum.def, aes(tree, loglik)) + 
  geom_line() + 
  facet_wrap(facets = ~panel.name) + 
  theme_classic() + 
  xlab('Iteration') + 
  ylab('Log likelihood')
pl.loglik.trace

ggarrange(plotlist = list(pl.nl.post, pl.loglik.trace), ncol = 1)



# set bornout and select iterations
remove.bornout <- function(df.sum, born.out, sample.per.chain, n.chain){
  df.split <- split(df.sum, df.sum$panel.name)
  idx.bornout.sm <- sort(unlist(outer(1:born.out, sample.per.chain*(0:(n.chain - 1)), '+')))
  df.split.red <- lapply(df.split, function(x) x[-idx.bornout.sm,] %>% 
                           mutate(new.idx = 1:(nrow(x) - length(idx.bornout.sm))))
  Reduce(rbind, df.split.red)
}




ggplot(df.sum.def, aes(x = nterm, y = loglik)) + 
  #geom_vline(xintercept = , color = 'grey', alpha = 0.3) + 
  geom_point() + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  theme_classic() + 
  ylab('log Likelihood') + 
  xlab('Number of terminal nodes')


ggplot(df.sum.def, aes(x = nterm, y = miss.rate)) + 
  #geom_vline(xintercept = , color = 'grey', alpha = 0.3) + 
  geom_point() + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  theme_classic() + 
  ylab('Missing rate') + 
  xlab('Number of terminal nodes')



ggplot(df.sum.def[order(df.sum.def$nterm),], 
       aes(x = miss.rate, y = loglik, color = nterm.cut)) + 
  #geom_vline(xintercept = , color = 'grey', alpha = 0.3) + 
  geom_point() + 
  labs(color = ~n[L]) + 
  scale_color_viridis_d() + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  theme_classic() + 
  ylab('Log Likelihood') + 
  xlab('Missing Rate')


ggplot(df.sum.def, aes(x = loglik, #y = after_stat(density), 
                       color = panel.name, fill = panel.name)) + 
  geom_histogram(bins = 50, alpha = 0.5) 


# best miss rate and loglik
vapply(df.split, \(x) min(x$miss.rate), 0)
vapply(df.split, \(x) max(x$loglik), 0)

#########################
## ALTERNATIVE MODELS ###
######################### 



model.list.alt <- list(mcmc_lb.medium.lb,
                       mcmc_lb.medium.hb,
                       mcmc_lb.small.lb,
                       mcmc_lb.small.hb)
names(model.list.alt) <- c('LB - o = 0.30, g = 0.5',
                           'LB - o = 0.30, g = 1.5',
                           'LB - o = 0.42, g = 0.5',
                           'LB - o = 0.42, g = 1.5')


# extract depth, number of terminal nodes, missing rate and loglik of all the trees
depth.df.alt <- apply_fun_models(fun_ = get_depth, 
                                 mcmc.list = model.list.alt,
                                 born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)
nterm.df.alt <- apply_fun_models(fun_ = get_num_terminal_nodes, 
                                 mcmc.list = model.list.alt, 
                                 born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)

rapid.miss.rate <- function(tree_top, X, Y){
  pp <- get_value_tree(tree_top, X)
  pred_ <- pp >= 0.5
  sum(pred_ != Y)
}

miss.rate.df.alt <- apply_fun_models(fun_ = \(x) rapid.miss.rate(tree_top = x, X = X_pred.norm, Y = Y_mal),
                                     mcmc.list = model.list.alt,
                                     born.out.pc = 250, 
                                     n.chain = n.chain_par, 
                                     sample.pc = n.iter_par)


like.df.alt <- apply_fun_models(fun_ = \(x) cart_log_lik_binary(tree_top = x, Y = Y_mal, 
                                                                X = X_pred.norm), 
                                mcmc.list = model.list.alt, 
                                born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)




df.sum.alt <- data.frame(tree = nterm.df.alt$x,
                         panel.name = nterm.df.alt$panel.name,
                         loglik = like.df.alt$y,
                         nterm = nterm.df.alt$y,
                         miss.rate = miss.rate.df.alt$y)


df.sum.alt$nterm.cut <- cut(df.sum.alt$nterm, breaks = c(3,7,10,13,25))


df.nl.prior.alt <- rbind(data.frame(x = 1:20, 
                                    y = prior.nterm(1:20, 0.3),
                                    panel.name = 'LB - o = 0.30, g = 0.5'),
                         data.frame(x = 1:20, 
                                    y = prior.nterm(1:20, 0.42),
                                    panel.name = 'LB - o = 0.42, g = 0.5'),
                         data.frame(x = 1:20, 
                                    y = prior.nterm(1:20, 0.3),
                                    panel.name = 'LB - o = 0.30, g = 1.5'),
                         data.frame(x = 1:20, 
                                    y = prior.nterm(1:20, 0.42),
                                    panel.name = 'LB - o = 0.42, g = 1.5'))

# create df for post quantities
df.split.alt <- split(df.sum.alt, df.sum.alt$panel.name)
post.nl.mean.alt <- vapply(df.split.alt, \(x) mean(x$nterm), 0)
post.nl.q0.025.alt <- vapply(df.split.alt, \(x) quantile(x$nterm, 0.025), 0)
post.nl.q0.975.alt <- vapply(df.split.alt, \(x) quantile(x$nterm, 0.975), 0)
df.post.nl.sum.alt <- data.frame(mean = post.nl.mean.alt, 
                                 low = post.nl.q0.025.alt, 
                                 up = post.nl.q0.975.alt,
                                 panel.name = names(df.split.alt))


factor.levs <- c( 'LB - o = 0.30, g = 0.5',
                  'LB - o = 0.42, g = 0.5',
                  'LB - o = 0.30, g = 1.5',
                  'LB - o = 0.42, g = 1.5')


factor.labs <- c( 'LB - o = 0.30, g = 0.5' = parse(text=latex2exp::TeX('$LB - \\omega = 0.3, \\gamma = 0.5$')),
                  'LB - o = 0.42, g = 0.5' = parse(text=latex2exp::TeX('$LB - \\omega = 0.42, \\gamma = 0.5$')),
                  'LB - o = 0.30, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.3, \\gamma = 1.5$')),
                  'LB - o = 0.42, g = 1.5'= parse(text=latex2exp::TeX('$LB - \\omega = 0.42, \\gamma = 1.5$')))


df.sum.alt$panel.name <- factor(df.sum.alt$panel.name,
                                levels = factor.levs,
                                labels = factor.labs,
                                ordered = TRUE)
df.nl.prior.alt$panel.name <- factor(df.nl.prior.alt$panel.name,
                                     levels = factor.levs,
                                     labels = factor.labs,
                                     ordered = TRUE)
df.post.nl.sum.alt$panel.name <- factor(df.post.nl.sum.alt$panel.name,
                                        levels = factor.levs,
                                        labels = factor.labs,
                                        ordered = TRUE)

pl.nl.post.alt <- ggplot(df.sum.alt, aes(x = nterm, y = after_stat(density))) + 
  geom_histogram(binwidth = 1, color = 'black', fill = 'white') + 
  geom_vline(xintercept = seq(1,21,by = 2), color = 'grey', alpha = 0.3) + 
  geom_line(data = df.nl.prior.alt, aes(x,y), colour = 'deepskyblue', linewidth = 0.4) +
  geom_point(data = df.nl.prior.alt, aes(x,y), colour = 'deepskyblue', size = 0.4) + 
  geom_vline(data = df.post.nl.sum.alt, aes(xintercept = mean)) +
  geom_vline(data = df.post.nl.sum.alt, aes(xintercept = low), linetype = 2) +
  geom_vline(data = df.post.nl.sum.alt, aes(xintercept = up), linetype = 2) +
  facet_wrap(facets = ~panel.name, ncol = 2, labeller = label_parsed) + 
  theme_classic() + 
  xlab('Number of terminal nodes') + 
  scale_x_continuous(breaks = seq(1,21,by = 2))
pl.nl.post.alt

ggplot(df.sum.alt, aes(x = nterm, y = loglik)) + 
  geom_point() + 
  facet_wrap(facets = ~panel.name, ncol = 2, labeller = label_parsed) + 
  theme_classic() + 
  xlab('Number of terminal nodes') +
  ylab('Log likelihood') 


ggplot(df.sum.alt, aes(x = nterm, y = miss.rate)) + 
  geom_point() + 
  facet_wrap(facets = ~panel.name, ncol = 2, labeller = label_parsed) + 
  theme_classic() + 
  xlab('Number of terminal nodes') +
  ylab('Missing rate') 

ggplot(df.sum.alt[order(df.sum.alt$nterm),], aes(x = miss.rate, y = loglik, color = nterm.cut)) + 
  geom_point() + 
  facet_wrap(facets = ~panel.name, ncol = 2, labeller = label_parsed) + 
  theme_classic() + 
  xlab('Missing rate') +
  ylab('Log likelihood') + 
  scale_color_viridis_d() + 
  labs(color = ~n[L])









