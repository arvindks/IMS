%% This is an real data analysis
%  Data: This is the Radiogenomics data from Arvind, 
%        which is a 533-by-48 matrix.
%  Purpose: 1. Compare cg with Hutchinson vs Quasi-Newton with Hutchinson
%              under different missing fraction
%
%           Comparison metrics inlucde:
%           a. BIC 
%           b. Runtime
%           c. MSE over missing entries, observed entries and all entries 
%
%  Date: 2016-05-03
%  Author: Liuyi Hu lhu@ncsu.edu
%
%% Load on the data
X = dlmread('S:\Documents\BMC\data\radgeno_big\Z.txt', ',', 1, 1);
% Laplacian matrix
row_laplacian = dlmread('S:\Documents\BMC\data\radgeno_big\Lr_15.txt',...
                        ',', 1, 1);
col_laplacian = dlmread('S:\Documents\BMC\data\radgeno_big\Lc_5.txt',...
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
result_runtime = zeros(n_repeat, 2*n_sz, n_miss_pct); 
result_bic = zeros(n_repeat, 2*n_sz, n_miss_pct);
result_mse_missing = zeros(n_repeat, 2*n_sz, n_miss_pct);
result_mse_obs = zeros(n_repeat, 2*n_sz, n_miss_pct);
result_mse_all = zeros(n_repeat, 2*n_sz, n_miss_pct);

% Loop
for miss_pct_idx = 1:n_miss_pct
    
    miss_pct = miss_pct_list(miss_pct_idx);
    fprintf('Current missing percentage: %f \n', miss_pct);
    
    for repeat_idx = 1:n_repeat
        
        fprintf('Current repeat_idx: %d \n',repeat_idx);
        rng(1991 + repeat_idx * 11);
        X_miss= missing_pattern(X,miss_pct, 'pattern', 'random');
        x = X(:);
        miss_idx = isnan(X_miss(:));
        
        for hut_sz_idx =1:n_sz
            sample_size = sample_size_list(hut_sz_idx);
            % qn with Hucthinson    
            tic;
            [X_complete,gamma_best,bic_best] = bmc_qn(X_miss,[0 0],...
                row_laplacian,col_laplacian,'approximate',true,...
                'sample_size',sample_size,'show', false);
            result_runtime(repeat_idx, 2*hut_sz_idx-1, miss_pct_idx) = toc;
            result_bic(repeat_idx, 2*hut_sz_idx-1, miss_pct_idx) = bic_best;
            z = X_complete(:);
            result_mse_missing(repeat_idx, 2*hut_sz_idx-1, miss_pct_idx) = ...
                norm(z(miss_idx) - x(miss_idx)).^2/sum(miss_idx);
            result_mse_obs(repeat_idx, 2*hut_sz_idx-1, miss_pct_idx) = ...
                norm(z(~miss_idx) - x(~miss_idx)).^2/sum(~miss_idx);
            result_mse_all(repeat_idx, 2*hut_sz_idx-1, miss_pct_idx) = ...
                norm(z - x).^2/length(miss_idx);
 
            % cg with Hutchinson
            tic;
            [X_complete,gamma_best,bic_best] = bmc_qn_cg(X_miss,[0 0],...
                row_laplacian,col_laplacian,'approximate',true,...
                'sample_size',sample_size,'show', false);
            result_runtime(repeat_idx, 2*hut_sz_idx,miss_pct_idx) = toc;
            result_bic(repeat_idx, 2*hut_sz_idx,miss_pct_idx) = bic_best;
            z = X_complete(:);
            result_mse_missing(repeat_idx, 2*hut_sz_idx, miss_pct_idx) = ...
                 norm(z(miss_idx) - x(miss_idx)).^2/sum(miss_idx);
            result_mse_obs(repeat_idx, 2*hut_sz_idx,miss_pct_idx) = ...
                  norm(z(~miss_idx) - x(~miss_idx)).^2/sum(~miss_idx);
            result_mse_all(repeat_idx, 2*hut_sz_idx, miss_pct_idx) = ...
                 norm(z - x).^2/length(miss_idx);
        end
    
    
    end
    
end

%% Save results
FileName = ['S:\Documents\BMC\results\','demo_radgeno2_ex1_',...
            datestr(now, 'mm-dd-HH-MM')];
FileName_runtime = [FileName,'_runtime.mat'];
FileName_BIC = [FileName,'_BIC.mat'];
FileName_MSE_missing = [FileName,'_MSE_missing.mat'];
FileName_MSE_obs = [FileName,'_MSE_obs.mat'];
FileName_MSE_all = [FileName,'_MSE_all.mat'];
save(FileName_runtime, 'result_runtime');
save(FileName_BIC, 'result_bic' );
save(FileName_MSE_missing, 'result_mse_missing' );
save(FileName_MSE_obs, 'result_mse_obs' );
save(FileName_MSE_all, 'result_mse_all' );

