%% This is an real data analysis
%  Data: This is the Mondrian data
%        which is a 370-by-380 matrix.
%  Purpose: Quasi-Neweton with Hutchinson of different sample size
%           using AIC as objective
%           Comparison metrics inlucdes:
%           a. AIC 
%           b. Runtime
%           c. MSE over missing entries, observed entries and all entries
%
%  Note: When calculating the MSE, divide the complete matrix and the
%        observed matrix values by 255
%  Date: 2016-09-19
%  
%  Author: Liuyi Hu lhu@ncsu.edu
%

clear;
close all;
%% Load on the data
filename = 'S:\Documents\BMC\data\';
A = imread([filename,'compositionabypietmondrian.jpg']);
B = rgb2gray(A);
X = double(B);
X_true = X;
% Laplacian matrix
row_laplacian = dlmread('S:\Documents\BMC\data\Mondrian\Lr_5.txt',...
                        ',', 1, 1);
col_laplacian = dlmread('S:\Documents\BMC\data\Mondrian\Lc_5.txt',...
                        ',', 1, 1);

%% Start experiments
% Input
loggamma0 = [0 0];
miss_pct_list = [0.1 0.3 0.5];
n_miss_pct = length(miss_pct_list);
n_repeat = 30;
sample_size_list = [5 10 50];
n_sz = length(sample_size_list); 
% Initialize to store results
result_runtime = zeros(n_repeat, n_sz, n_miss_pct); 
result_aic = zeros(n_repeat, n_sz, n_miss_pct);
result_mse_missing = zeros(n_repeat, n_sz, n_miss_pct);
result_mse_obs = zeros(n_repeat, n_sz, n_miss_pct);
result_mse_all = zeros(n_repeat, n_sz, n_miss_pct);

sigma = 0.5*255;
% Loop
for miss_pct_idx = 1:n_miss_pct
    
    miss_pct = miss_pct_list(miss_pct_idx);
    fprintf('Current missing percentage: %f \n', miss_pct);
    
    for repeat_idx = 1:n_repeat
        
        fprintf('Current repeat_idx: %d \n',repeat_idx);
        rng(1991 + repeat_idx * 11);
        % Add noise
        X = X_true + sigma * randn(size(X_true));
        % Add missing
        X_miss= missing_pattern(X, miss_pct, 'pattern', 'random');
        x = X(:);
        miss_idx = isnan(X_miss(:));
        for hut_sz_idx =1:n_sz
            sample_size = sample_size_list(hut_sz_idx);
            % Quasi-Newton direct with Hucthinson    
            tic;
            [X_complete,gamma_best,aic_best] = bmc_qn(X_miss,[0 0],...
                row_laplacian,col_laplacian,'approximate',true,...
                'sample_size',sample_size,'show', false,'criterion', 'AIC');
            result_runtime(repeat_idx, hut_sz_idx, miss_pct_idx) = toc;
            result_aic(repeat_idx, hut_sz_idx, miss_pct_idx) = aic_best;
            z = X_complete(:);
            result_mse_missing(repeat_idx, hut_sz_idx, miss_pct_idx) = ...
                norm(z(miss_idx) - x(miss_idx)).^2/sum(miss_idx) / 255^2;
            result_mse_obs(repeat_idx, hut_sz_idx, miss_pct_idx) = ...
                norm(z(~miss_idx) - x(~miss_idx)).^2/sum(~miss_idx) / 255^2;
            result_mse_all(repeat_idx, hut_sz_idx, miss_pct_idx) = ...
                norm(z - x).^2/length(miss_idx) / 255^2;
 
        end

    
    
    end
    
end

%% Save results
FileName = ['S:\Documents\BMC\results\','demo_Mondrian_ex2_AIC_',...
            datestr(now, 'mm-dd-HH-MM')];
FileName_runtime = [FileName,'_runtime.mat'];
FileName_AIC = [FileName,'_AIC.mat'];
FileName_MSE_missing = [FileName,'_MSE_missing.mat'];
FileName_MSE_obs = [FileName,'_MSE_obs.mat'];
FileName_MSE_all = [FileName,'_MSE_all.mat'];
save(FileName_runtime, 'result_runtime');
save(FileName_AIC, 'result_aic' );
save(FileName_MSE_missing, 'result_mse_missing' );
save(FileName_MSE_obs, 'result_mse_obs' );
save(FileName_MSE_all, 'result_mse_all' );






