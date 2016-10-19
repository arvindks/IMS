%% This script saves data for the plot of "Figure5_experiments.m" to .csv format
%
%  Date: 2016-05-27, 2016-09-14
%  Author: Eric Chi, eric_chi@ncsu.edu
%
%% Load results BIC
% Edit the next line for local path info.
prefix = '/Users/ecchi/Dropbox/Work/Research/00_Active/00_NCSU/01_BMC/results/';
%suffix = 'Figure5_experiments_05-25-20-13';
suffix = 'Figure5_experiments_08-25-20-54';
FileName = [prefix, suffix];
FileName_runtime = [FileName,'_runtime.mat'];
FileName_BIC = [FileName,'_BIC.mat'];
FileName_MSE_missing = [FileName,'_MSE_missing.mat'];
FileName_MSE_obs = [FileName,'_MSE_obs.mat'];
FileName_MSE_all = [FileName,'_MSE_all.mat'];
load(FileName_runtime, 'result_runtime');
load(FileName_BIC, 'result_bic' );
load(FileName_MSE_missing, 'result_mse_missing' );
load(FileName_MSE_obs, 'result_mse_obs' );
load(FileName_MSE_all, 'result_mse_all' );

%% Input
n_finer = 3;
miss_pct_list = [0.1 0.3 0.5];
n_miss_pct = length(miss_pct_list);
n_repeat = 30;

%% Filename Out
% Save information on the fraction of missing entries
FileNameOut_setup = [FileName,'_setup_fraction.csv'];
csvwrite(FileNameOut_setup, miss_pct_list);

% Save Run Times
for i = 1:n_miss_pct
   FileNameOut_runtime = [FileName,'_runtime_',sprintf('%g',miss_pct_list(i)) ,'.csv'];
   csvwrite(FileNameOut_runtime, result_runtime(:,:,i));
end

% Save BIC scores
for i = 1:n_miss_pct
   FileNameOut_bic = [FileName,'_bic_',sprintf('%g',miss_pct_list(i)) ,'.csv'];
   csvwrite(FileNameOut_bic, result_bic(:,:,i));
end

% Save MSE over missing entries
for i = 1:n_miss_pct
   FileNameOut_mse_missing = [FileName,'_mse_missing_',sprintf('%g',miss_pct_list(i)) ,'.csv'];
   csvwrite(FileNameOut_mse_missing, result_mse_missing(:,:,i));
end

% Save MSE over observed entries
for i = 1:n_miss_pct
   FileNameOut_mse_obs = [FileName,'_mse_obs_',sprintf('%g',miss_pct_list(i)) ,'.csv'];
   csvwrite(FileNameOut_mse_obs, result_mse_obs(:,:,i));
end

%% Load results AIC
prefix = '/Users/ecchi/Dropbox/Work/Research/00_Active/00_NCSU/01_BMC/results/';
suffix = 'numerical_ex6_AIC_08-28-02-07';
FileName = [prefix, suffix];
FileName_runtime = [FileName,'_runtime.mat'];
FileName_AIC = [FileName,'_AIC.mat'];
FileName_MSE_missing = [FileName,'_MSE_missing.mat'];
FileName_MSE_obs = [FileName,'_MSE_obs.mat'];
FileName_MSE_all = [FileName,'_MSE_all.mat'];
load(FileName_runtime, 'result_runtime');
load(FileName_AIC, 'result_aic' );
load(FileName_MSE_missing, 'result_mse_missing' );
load(FileName_MSE_obs, 'result_mse_obs' );
load(FileName_MSE_all, 'result_mse_all' );

%% Input
n_finer = 3;
miss_pct_list = [0.1 0.3 0.5];
n_miss_pct = length(miss_pct_list);
n_repeat = 30;

%% Filename Out
% Save information on the fraction of missing entries
FileNameOut_setup = [FileName,'_setup_fraction.csv'];
csvwrite(FileNameOut_setup, miss_pct_list);

% Save Run Times
for i = 1:n_miss_pct
   FileNameOut_runtime = [FileName,'_runtime_',sprintf('%g',miss_pct_list(i)) ,'.csv'];
   csvwrite(FileNameOut_runtime, result_runtime(:,:,i));
end

% Save AIC scores
for i = 1:n_miss_pct
   FileNameOut_aic = [FileName,'_aic_',sprintf('%g',miss_pct_list(i)) ,'.csv'];
   csvwrite(FileNameOut_aic, result_aic(:,:,i));
end

% Save MSE over missing entries
for i = 1:n_miss_pct
   FileNameOut_mse_missing = [FileName,'_mse_missing_',sprintf('%g',miss_pct_list(i)) ,'.csv'];
   csvwrite(FileNameOut_mse_missing, result_mse_missing(:,:,i));
end

% Save MSE over observed entries
for i = 1:n_miss_pct
   FileNameOut_mse_obs = [FileName,'_mse_obs_',sprintf('%g',miss_pct_list(i)) ,'.csv'];
   csvwrite(FileNameOut_mse_obs, result_mse_obs(:,:,i));
end