library(ggplot2)
library(ggthemes)
library(reshape2)

## Manually adjust
#filename <- 'Figure3_experiments_05-24-22-03'
filename <- 'Figure3_experiments_08-24-19-42'
method_order <- c(6,1:5,7)

sample_size <- read.table(paste0(filename,'_setup_sample_size.csv'),sep=',',header=FALSE)
#sample_size <- paste('Sample Size =', sample_size)
sample_size <- paste0('QN (H', sample_size,")")
n_sample_size <- length(sample_size)

missing_fraction <- read.table(paste0(filename,'_setup_fraction.csv'),sep=',',header=FALSE)
n_missing_fraction <- length(missing_fraction)

mse_obs <- mse_missing <- bic <- runtime <- matrix(NA, 0, n_sample_size+2)

for (i in 1:n_missing_fraction) {
  d_runtime <- read.table(paste0(filename,'_runtime_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_runtime[,n_sample_size+2] <- as.double(missing_fraction[i])
  d_runtime[,n_sample_size+2] <- as.factor(d_runtime[,n_sample_size+2])
  runtime <- rbind(runtime, d_runtime)
  
  d_bic <- read.table(paste0(filename,'_bic_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_bic[,n_sample_size+2] <- as.double(missing_fraction[i])
  d_bic[,n_sample_size+2] <- as.factor(d_bic[,n_sample_size+2])
  bic <- rbind(bic, d_bic)
  
  d_mse_missing <- read.table(paste0(filename,'_mse_missing_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_mse_missing[,n_sample_size+2] <- as.double(missing_fraction[i])
  d_mse_missing[,n_sample_size+2] <- as.factor(d_mse_missing[,n_sample_size+2])
  mse_missing <- rbind(mse_missing, d_mse_missing)
  
  d_mse_obs <- read.table(paste0(filename,'_mse_obs_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_mse_obs[,n_sample_size+2] <- as.double(missing_fraction[i])
  d_mse_obs[,n_sample_size+2] <- as.factor(d_mse_obs[,n_sample_size+2])
  mse_obs <- rbind(mse_obs, d_mse_obs)
}

names(mse_obs) <- names(mse_missing) <- names(bic) <- names(runtime) <- c(sample_size, 'QN (E)', 'Fraction Missing')

runtime <- runtime[,method_order]
runtime <- melt(runtime, id.vars='Fraction Missing')
runtime$value <- log10(runtime$value)
#q <- ggplot(data=runtime, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) + scale_y_log10()

bic <- bic[,method_order]
bic <- melt(bic, id.vars='Fraction Missing')
n <- 1e2
bic$value <- bic$value / (n**2)
#bic$value <- log10(bic$value)
#q <- ggplot(data=bic, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) #+ scale_y_log10()

mse_missing <- mse_missing[,method_order]
mse_missing <- melt(mse_missing, id.vars='Fraction Missing')
#mse_missing$value <- log10(mse_missing$value)
#q <- ggplot(data=mse_missing, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) #+ scale_y_log10()

mse_obs <- mse_obs[,method_order]
mse_obs <- melt(mse_obs, id.vars='Fraction Missing')
#mse_obs$value <- log10(mse_obs$value)
#q <- ggplot(data=mse_obs, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) #+ scale_y_log10()
## Create a single figures
runtime$metric <- 'log10(Time in sec)'
bic$metric <- 'BIC'
mse_missing$metric <- 'MSE over missing entries'
mse_obs$metric <- 'MSE over observed entries'

df <- rbind(mse_obs,rbind(mse_missing,rbind(bic, runtime)))
df$metric <- ordered(df$metric, levels = c('log10(Time in sec)', 'BIC', 'MSE over missing entries', 'MSE over observed entries'))

q <- ggplot(data=df, aes(x=`Fraction Missing`, y=value))
q <- q + geom_boxplot(aes(fill=variable)) + scale_fill_brewer(guide = guide_legend(title=NULL),palette="RdYlBu")
q + facet_wrap(~metric, scales='free_y') + ylab('') + theme(legend.text = element_text(size = 13), strip.text = element_text(size = 14), axis.text = element_text(size=14), axis.title=element_text(size=14))

golden_ratio <- 1.61803398875
height <- 7
filename <- 'numerical_ex4.pdf'
ggsave(filename, height=height, width=golden_ratio*height)


## Part II: AIC

## Manually adjust
rm(list=ls())
filename <- 'Figure3_experiments_AIC_08-26-22-15'
method_order <- c(6,1:5,7)

sample_size <- read.table(paste0(filename,'_setup_sample_size.csv'),sep=',',header=FALSE)
#sample_size <- paste('Sample Size =', sample_size)
sample_size <- paste0('QN (H', sample_size,")")
n_sample_size <- length(sample_size)

missing_fraction <- read.table(paste0(filename,'_setup_fraction.csv'),sep=',',header=FALSE)
n_missing_fraction <- length(missing_fraction)

mse_obs <- mse_missing <- aic <- runtime <- matrix(NA, 0, n_sample_size+2)

for (i in 1:n_missing_fraction) {
  d_runtime <- read.table(paste0(filename,'_runtime_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_runtime[,n_sample_size+2] <- as.double(missing_fraction[i])
  d_runtime[,n_sample_size+2] <- as.factor(d_runtime[,n_sample_size+2])
  runtime <- rbind(runtime, d_runtime)
  
  d_aic <- read.table(paste0(filename,'_aic_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_aic[,n_sample_size+2] <- as.double(missing_fraction[i])
  d_aic[,n_sample_size+2] <- as.factor(d_aic[,n_sample_size+2])
  aic <- rbind(aic, d_aic)
  
  d_mse_missing <- read.table(paste0(filename,'_mse_missing_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_mse_missing[,n_sample_size+2] <- as.double(missing_fraction[i])
  d_mse_missing[,n_sample_size+2] <- as.factor(d_mse_missing[,n_sample_size+2])
  mse_missing <- rbind(mse_missing, d_mse_missing)
  
  d_mse_obs <- read.table(paste0(filename,'_mse_obs_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_mse_obs[,n_sample_size+2] <- as.double(missing_fraction[i])
  d_mse_obs[,n_sample_size+2] <- as.factor(d_mse_obs[,n_sample_size+2])
  mse_obs <- rbind(mse_obs, d_mse_obs)
}

names(mse_obs) <- names(mse_missing) <- names(aic) <- names(runtime) <- c(sample_size, 'QN (E)', 'Fraction Missing')

runtime <- runtime[,method_order]
runtime <- melt(runtime, id.vars='Fraction Missing')
runtime$value <- log10(runtime$value)
#q <- ggplot(data=runtime, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) + scale_y_log10()

aic <- aic[,method_order]
aic <- melt(aic, id.vars='Fraction Missing')
n <- 1e2
aic$value <- aic$value / (n**2)
#bic$value <- log10(bic$value)
#q <- ggplot(data=bic, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) #+ scale_y_log10()

mse_missing <- mse_missing[,method_order]
mse_missing <- melt(mse_missing, id.vars='Fraction Missing')
#mse_missing$value <- log10(mse_missing$value)
#q <- ggplot(data=mse_missing, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) #+ scale_y_log10()

mse_obs <- mse_obs[,method_order]
mse_obs <- melt(mse_obs, id.vars='Fraction Missing')
#mse_obs$value <- log10(mse_obs$value)
#q <- ggplot(data=mse_obs, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) #+ scale_y_log10()
## Create a single figures
runtime$metric <- 'log10(Time in sec)'
aic$metric <- 'AIC'
mse_missing$metric <- 'MSE over missing entries'
mse_obs$metric <- 'MSE over observed entries'

df <- rbind(mse_obs,rbind(mse_missing,rbind(aic, runtime)))
df$metric <- ordered(df$metric, levels = c('log10(Time in sec)', 'AIC', 'MSE over missing entries', 'MSE over observed entries'))

q <- ggplot(data=df, aes(x=`Fraction Missing`, y=value))
q <- q + geom_boxplot(aes(fill=variable)) + scale_fill_brewer(guide = guide_legend(title=NULL),palette="RdYlBu")
q + facet_wrap(~metric, scales='free_y') + ylab('') + theme(legend.text = element_text(size = 13), strip.text = element_text(size = 14), axis.text = element_text(size=14), axis.title=element_text(size=14))

golden_ratio <- 1.61803398875
height <- 7
filename <- 'numerical_ex4_AIC.pdf'
ggsave(filename, height=height, width=golden_ratio*height)
