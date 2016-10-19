%% This is a numerical experiment
%  Purpose: compare performance using different tuning methods
%           1. Fminunc with exact solve
%           2. AIC grid search with exact solve
%
%           Comparison metrics inlucdes:
%           a. AIC 
%           b. Runtime
%           c. MSE over missing entries, observed entries and all entries
%
%  Date: 2016-08-26
%  Author: Liuyi Hu lhu@ncsu.edu

%% Set up the missing matrix 
blocks = [25 25 25 25];
%blocks = [5 5 5 5];
X11 = -25 * ones(blocks(1), blocks(1)); %1st block
X12 = 10 * ones(blocks(2), blocks(2));  %2nd blcok
X21 = -10 * ones(blocks(3), blocks(3)); %3rd block
X22 = 25 * ones(blocks(4), blocks(4));  %4th block
X = [X11 X12; X21 X22];

% shuffle the matrix
[n, p] = size(X);
% add some noise
X_noise = X + normrnd(0,1,n,p);
% Laplacian
row_laplacian = laplacian(n, [blocks(1), blocks(2)]);
col_laplacian = laplacian(p, [blocks(3), blocks(4)]);

%% Start experiments
% Input
init_grid_n = 10; % initial grid points
gamma0 = [1 1];
loggamma0 = [0 0];
n_finer = 3;
miss_pct_list = [0.1 0.3 0.5];
n_miss_pct = length(miss_pct_list);
n_repeat = 30;

% Initialize to store results
result_runtime = zeros(n_repeat, n_finer+2, n_miss_pct); 
result_aic = zeros(n_repeat, n_finer+2, n_miss_pct);
result_mse_missing = zeros(n_repeat, n_finer+2, n_miss_pct);
result_mse_obs = zeros(n_repeat, n_finer+2, n_miss_pct);
result_mse_all = zeros(n_repeat, n_finer+2, n_miss_pct);

% Loop
for miss_pct_idx = 1:n_miss_pct
    
    miss_pct = miss_pct_list(miss_pct_idx);
    fprintf('Current missing percentage: %f \n', miss_pct);
    
    for repeat_idx = 1:n_repeat
        
        fprintf('Current repeat_idx: %d \n',repeat_idx);
        rng(1991 + repeat_idx * 11);
        X_miss= missing_pattern(X_noise, miss_pct, 'pattern', 'random');
        %x = X_noise(:);
        x = X(:);
        miss_idx = isnan(X_miss(:));
        
        for finer_idx =1:n_finer
            
            row_gamma_list = exp(linspace(-9, 1,... 
                   init_grid_n+(init_grid_n-1)*(2^(finer_idx-1)-1)));
            col_gamma_list = exp(linspace(-9, 1,...
                   init_grid_n+(init_grid_n-1)*(2^(finer_idx-1)-1)));
            tic;
            [X_complete,~,aic_best,aic_grid] = bmc_grid(X_miss, ...
                                      row_gamma_list,col_gamma_list,...
                                      row_laplacian,col_laplacian,...
                                      'criterion', 'AIC');
                                  
            result_runtime(repeat_idx, finer_idx,miss_pct_idx) = toc;
            result_aic(repeat_idx, finer_idx,miss_pct_idx) = aic_best;
            
            z = X_complete(:);
            
            result_mse_missing(repeat_idx,finer_idx,miss_pct_idx) = ...
                norm(z(miss_idx) - x(miss_idx)).^2/sum(miss_idx);
            result_mse_obs(repeat_idx, finer_idx,miss_pct_idx) = ...
                norm(z(~miss_idx) - x(~miss_idx)).^2/sum(~miss_idx);
            result_mse_all(repeat_idx,finer_idx,miss_pct_idx) = ...
                norm(z - x).^2/length(miss_idx);
        
        end
        
        % fminunc
        tic;
        [X_complete,gamma_best,aic_best] = bmc_qn(X_miss,[0 0],...
                row_laplacian,col_laplacian, 'criterion', 'AIC');
        result_runtime(repeat_idx, n_finer+1,miss_pct_idx) = toc;
        result_aic(repeat_idx, n_finer+1,miss_pct_idx) = aic_best;
        z = X_complete(:);
        result_mse_missing(repeat_idx, n_finer+1, miss_pct_idx) = ...
                norm(z(miss_idx) - x(miss_idx)).^2/sum(miss_idx);
        result_mse_obs(repeat_idx, n_finer+1,miss_pct_idx) = ...
                norm(z(~miss_idx) - x(~miss_idx)).^2/sum(~miss_idx);
        result_mse_all(repeat_idx, n_finer+1, miss_pct_idx) = ...
                norm(z - x).^2/length(miss_idx);
            
         % fminunc with Hutchinson, sample size = 5
        tic;
        [X_complete,~,aic_best] = bmc_qn(X_miss,[0 0],...
                row_laplacian,col_laplacian,'approximate',true,...
                'sample_size',5, 'criterion', 'AIC');
        result_runtime(repeat_idx, n_finer+2,miss_pct_idx) = toc;
        result_aic(repeat_idx, n_finer+2,miss_pct_idx) = aic_best;
        z = X_complete(:);
        result_mse_missing(repeat_idx, n_finer+2, miss_pct_idx) = ...
                norm(z(miss_idx) - x(miss_idx)).^2/sum(miss_idx);
        result_mse_obs(repeat_idx, n_finer+2,miss_pct_idx) = ...
                norm(z(~miss_idx) - x(~miss_idx)).^2/sum(~miss_idx);
        result_mse_all(repeat_idx, n_finer+2, miss_pct_idx) = ...
                norm(z - x).^2/length(miss_idx);       
    
    end
    
end

%% Save results
FileName = ['S:\Documents\BMC\results\','numerical_ex6_AIC_',...
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


