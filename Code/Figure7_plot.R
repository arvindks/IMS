library(ggplot2)
library(ggthemes)
library(reshape2)

## Manually adjust
#filename <- 'demo_radgeno2_ex1_05-10-22-35'
filename <- 'demo_radgeno2_ex1_08-26-15-21'
#method_order <- c(1,3,5,2,4,6,7)

sample_size <- read.table(paste0(filename,'_setup_sample_size.csv'),sep=',',header=FALSE)
#sample_size <- paste('Sample Size =', sample_size)
sample_size <- c(paste0('QN (H', sample_size,")"), paste0('CG (H', sample_size,")"))
n_sample_size <- length(sample_size)

missing_fraction <- read.table(paste0(filename,'_setup_fraction.csv'),sep=',',header=FALSE)
n_missing_fraction <- length(missing_fraction)

mse_obs <- mse_missing <- bic <- runtime <- matrix(NA, 0, 2*n_sample_size)

for (i in 1:n_missing_fraction) {
  d_runtime <- read.table(paste0(filename,'_runtime_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_runtime <- d_runtime[,c(1,3,5,2,4,6)]
  d_runtime[,n_sample_size+1] <- as.double(missing_fraction[i])
  d_runtime[,n_sample_size+1] <- as.factor(d_runtime[,n_sample_size+1])
  runtime <- rbind(runtime, d_runtime)
  
  d_bic <- read.table(paste0(filename,'_bic_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_bic <- d_bic[,c(1,3,5,2,4,6)]
  d_bic[,n_sample_size+1] <- as.double(missing_fraction[i])
  d_bic[,n_sample_size+1] <- as.factor(d_bic[,n_sample_size+1])
  bic <- rbind(bic, d_bic)
  
  d_mse_missing <- read.table(paste0(filename,'_mse_missing_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_mse_missing <- d_mse_missing[,c(1,3,5,2,4,6)]
  d_mse_missing[,n_sample_size+1] <- as.double(missing_fraction[i])
  d_mse_missing[,n_sample_size+1] <- as.factor(d_mse_missing[,n_sample_size+1])
  mse_missing <- rbind(mse_missing, d_mse_missing)
  
  d_mse_obs <- read.table(paste0(filename,'_mse_obs_',missing_fraction[i],'.csv'),sep=',',header=FALSE)
  d_mse_obs <- d_mse_obs[,c(1,3,5,2,4,6)]
  d_mse_obs[,n_sample_size+1] <- as.double(missing_fraction[i])
  d_mse_obs[,n_sample_size+1] <- as.factor(d_mse_obs[,n_sample_size+1])
  mse_obs <- rbind(mse_obs, d_mse_obs)
}

names(mse_obs) <- names(mse_missing) <- names(bic) <- names(runtime) <- c(sample_size, 'Fraction Missing')

#runtime <- runtime[,method_order]
runtime <- melt(runtime, id.vars='Fraction Missing')
runtime$value <- log10(runtime$value)
#q <- ggplot(data=runtime, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) + scale_y_log10()

#bic <- bic[,method_order]
bic <- melt(bic, id.vars='Fraction Missing')
n <- 1e2
bic$value <- bic$value / (n**2)
#bic$value <- log10(bic$value)
#q <- ggplot(data=bic, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) #+ scale_y_log10()

#mse_missing <- mse_missing[,method_order]
mse_missing <- melt(mse_missing, id.vars='Fraction Missing')
#mse_missing$value <- log10(mse_missing$value)
#q <- ggplot(data=mse_missing, aes(x=`Fraction Missing`, y=value))
#q + geom_boxplot(aes(fill=variable)) #+ scale_y_log10()

#mse_obs <- mse_obs[,method_order]
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
filename <- 'radgeno.pdf'
ggsave(filename, height=height, width=golden_ratio*height)