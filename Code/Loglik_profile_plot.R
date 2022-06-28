require(ggplot2)
require(cowplot)
require(stats)
require(gtools)
# Load results of log-likelihood profile analyses
load(file="Liklihood_profile_R.rdata")
df_Var_loglik_R = df_VarPD_loglik
load(file="Liklihood_profile_S.rdata")
df_Var_loglik_S = df_VarPD_loglik
load(file="Liklihood_profile_P.rdata")
df_Var_loglik_P = df_VarPD_loglik
#
df_loglik = data.frame(Parameter = c(rep("Recruitment",11),rep("Survival",11),rep("Detection",11)),
                       Prctn_Var = c(df_Var_loglik_R$VarPD_ppn*100,
                                     df_Var_loglik_S$VarPD_ppn*100,
                                     df_Var_loglik_P$VarPD_ppn*100),
                       Loglik = c(df_Var_loglik_R$elpd_loo,
                                  df_Var_loglik_S$elpd_loo,
                                  df_Var_loglik_P$elpd_loo))
df_loglik$Parameter = factor(df_loglik$Parameter)
#
ggplot(df_loglik,aes(x=Prctn_Var,y=Loglik,group=Parameter,fill=Parameter,color=Parameter)) +
  # geom_point() +
  geom_smooth() +
  labs(x="Percent of estimated sigma",y="Expected log predictive density (elpd)") +
  theme_classic()
