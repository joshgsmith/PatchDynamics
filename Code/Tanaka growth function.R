# setwd('/Users/Joshua/Documents/Research/Papers/chapter 1 - Urchin behavioral shift/Analyses/H1 - Size frequency analysis/Process Model/My model/Final Model')
#required packages
library(vcdExtra)
library(reshape2)
library(ggplot2)
library(dplyr)

DatU.raw <- read.csv('./Data/Size_fq.csv') #Urchin size frequency
# DatS <- read.csv('./Data/Stat_sum.csv') #Mean density by size class -- Stat_sum.csv

# size_fq <- read.csv("./Data/size_year.csv")
# DatU.raw <- size_fq %>%
#   dplyr::group_by(size) %>%
#   dplyr::summarize(Freq = sum(count))
#attach(DatU.raw)
#attach(DatS)
#colnames(DatS)[1] = "Year"
colnames(DatU.raw)[1] = "UrchSize"
FrequencyTable <- DatU.raw
DatU <- expand.dft(FrequencyTable, freq="Freq")

# From Ebert, T. A. (2010). Demographic patterns of the purple sea urchin Strongylocentrotus purpuratus along a latitudinal gradient, 1985–1987. Marine Ecology Progress Series, 406, 105-120. - use Bodega Bay for growth params
# Paramaters for Tanaka function:
# 
# NOTE: used figure 7 from Ebert to extract 95% CI for growth increment at 3 cm for Bodega (0.05cm)
# with sample size of 141, back calculated SD based on 95% CI = +/- 1.96 * SE/sqrt(n)
# Sig = .05*sqrt(141)/1.96
# gives SD of growth increment of .575 cm as 0.3 cm (~ CV of 0.3/.525 = 0.57)

a =  0.558
d =  1.432
f =  1.3
CV = .57
Sig = 0.3
reps = 500
# Define size classes:
Size = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) 
E = exp(sqrt(f)*((Size-.2)-d))
Size2 = 1/sqrt(f)*log(2*f*(E/(4*f) - a/E +1) + 2*sqrt((f^2)*(E/(4*f) - a/E +1)^2 + f*a )) + d 
Incsize = pmax(.1,Size2 - (Size-.2))
plot(Size,Incsize)
Dt_det = DatU$UrchSize[DatU$UrchSize<10] 
Inc_det = numeric(length = length(Dt_det))
for(i in 1:9){
  ii = which(Dt_det==Size[i])
  Inc_det[ii] = Incsize[i]
}
b = Inc_det / Sig^2
Gr_rnd = numeric(length = length(Dt_det))
for (r in 1:reps){
  Dt = Dt_det + runif(length(Dt_det), -.4999, .4999)
  Dt1 = round(Dt + rgamma(length(Dt),Inc_det*b,b))
  Gr_rnd = Gr_rnd + pmax(0,pmin(1,Dt1 - Dt_det))
}
df = data.frame(Sz = factor(Dt_det), GrP = Gr_rnd/reps)
G_hat = df %>% dplyr::group_by(Sz) %>% dplyr::summarise(mean = mean(GrP),
                                                        sd = sd(GrP))
print('G prob values: ')
noquote(paste(format(G_hat$mean,digits=3),collapse = ','))

# Box plot of growth prob
plt_Gr1 = ggplot(df,aes(x=Sz,y=GrP,group=Sz)) +
  geom_boxplot() +
  labs(x="Size class",y="Growth transition probability") +
  theme_classic()
print(plt_Gr1)

# Plot stochastic urchin growth increments scatter plot 
Size = seq(.1,10,by=.1) 
E = exp(sqrt(f)*((Size)-d))
Size2 = 1/sqrt(f)*log(2*f*(E/(4*f) - a/E +1) + 2*sqrt((f^2)*(E/(4*f) - a/E +1)^2 + f*a )) + d 
Incsize = pmax(.1,Size2 - (Size))
smth = smooth.spline(Size,Incsize) 
Dt = Dt_det + runif(length(Dt_det), -.4999, .4999)
E = exp(sqrt(f)*(Dt-d))
Dt1_det = 1/sqrt(f)*log(2*f*(E/(4*f) - a/E +1) + 2*sqrt((f^2)*(E/(4*f) - a/E +1)^2 + f*a )) + d
Inc_det = pmax(.2,Dt1_det - Dt)
b = Inc_det / Sig^2
Inc_stoch = rgamma(length(Dt),Inc_det*b,b)
ii = which(Dt_det == 1 | Dt_det>6)
ii = c(ii, sample(which(Dt_det == 2 | Dt_det == 6),1000))
ii = c(ii, sample(which(Dt_det > 2 | Dt_det<6),3000))
df2 = data.frame(Size=Dt[ii],Growth_increment=Inc_stoch[ii])
plt_Gr2 = ggplot(df2,aes(x=Size,y=Growth_increment)) +
  geom_point(shape=16,color="darkgrey") +
  geom_line(data=data.frame(Size=Size,Gr_mean=smth$y),aes(x=Size,y=Gr_mean),
            size=1.1,color="blue") +
  labs(x="Size (cm)",y="Growth increment (cm)") +
  theme_classic()
print(plt_Gr2)

# These are the outputs from the growth model, stage-specific transition probabilities
growth_matrix <- matrix(
  c(
    1-G_hat[1], 0, 0, 0, 0, 0, 0, 0, 0, 0,
    G_hat[1], 1-G_hat[2], 0, 0, 0, 0, 0, 0, 0, 0,
    0, G_hat[2], 1-G_hat[3], 0, 0, 0, 0, 0, 0, 0,
    0, 0, G_hat[3], 1-G_hat[4], 0, 0, 0, 0, 0, 0,
    0, 0, 0, G_hat[4], 1-G_hat[5], 0 ,0 ,0 ,0 ,0,
    0, 0, 0, 0, G_hat[5], 1-G_hat[6], 0, 0, 0, 0,
    0, 0, 0, 0, 0, G_hat[6], 1-G_hat[7],0 ,0, 0,
    0, 0, 0, 0, 0, 0, G_hat[7], 1-G_hat[8], 0, 0,
    0, 0, 0, 0, 0, 0, 0, G_hat[8], 1-G_hat[9], 0,
    0, 0, 0, 0, 0, 0, 0, 0, G_hat[9], 0),
  nrow=10, ncol=10, byrow=T
)


