########################################################
## Library
########################################################
library(rstan)
library(tidyverse)
library(plyr)
library(coda)
library(bayesplot)
library(survival)
library(purrr)
library(survminer)
library(tidyr)
library(broom)
########################################################
## Functions
########################################################
MCMC_traceplot <- function(fit, para_names = NULL, nb_chains = 3, burnin = TRUE){
  if(burnin ==TRUE){
    nb_para <- length(para_names)
    #print(nb_para)
    names(fit)[1:nb_para] <- para_names # you can also add  "loglik", "lp__" as other "parameters"
    p1 <- bayesplot::mcmc_trace(rstan::extract(fit, inc_warmup = TRUE, permuted = FALSE)[,1:nb_chains,1:nb_para],
                                n_warmup = 1000,
                                facet_args = list(nrow =nb_para, labeller = label_parsed))
    
    p2 <- p1 + bayesplot::facet_text(size = 15)+
      theme_bw()+
      theme(legend.position = "bottom")+
      scale_x_continuous(expand = c(0, 10))+
      scale_color_manual(values=c("#FF7F00", "#CB92E8","#5B1584")[1:nb_chains])
  }else{
    mcmc <- as.array(fit)[,,1:nb_para]
    p1 <- bayesplot::mcmc_trace(mcmc,
                     n_warmup = 0,
                     facet_args = list(nrow = 2, labeller = label_parsed))
    
    p2 <- p1 + bayesplot::facet_text(size = 15)+
      theme_bw()+
      scale_x_continuous(breaks = seq(0, 1000, 200), labels = seq(1000, 2000, 200))+
      scale_color_manual(values=c("#FF7F00", "#CB92E8","#5B1584" ))
  }

  return(p2)
}
