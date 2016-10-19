%% This is a demo of solving bmc using different methods

%% raw X 
blocks = [5 5 5 5];
X11 = -25 * ones(blocks(1), blocks(1)); %1st block
X12 = 10 * ones(blocks(2), blocks(2));  %2nd blcok
X21 = -10 * ones(blocks(3), blocks(3)); %3rd block
X22 = 25 * ones(blocks(4), blocks(4));  %4th block
X = [X11 X12; X21 X22];
%add noise
[n, p] = size(X);
X_noise = X + normrnd(0,1,n,p);
% impose missing
X_miss_random = missing_pattern(X_noise, 0.2, 'pattern','random');
% plot
imagesc(X_miss_random);
% Laplacian matrix
row_laplacian = laplacian(n, [blocks(1), blocks(2)]);
col_laplacian = laplacian(p, [blocks(3), blocks(4)]);

%% Nelder-Mead with exact solve
% BIC as the objective
[X_complete,gamma_best,bic_best] = bmc_nm(X_miss_random,[1 1],...
row_laplacian,col_laplacian,'approximate',false,'show',false);
%imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum BIC is %d .\n',bic_best);
% AIC as the objective
[X_complete,gamma_best,aic_best] = bmc_nm(X_miss_random,[1 1],...
row_laplacian,col_laplacian,'approximate',false,'show',false, ...
'criterion', 'AIC');
%imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum AIC is %d .\n',aic_best);

%% Nelder-Mead with Hutchinson estimator
% NOTE: This sometimes cannot converge
% BIC as the objective
[X_complete,gamma_best,bic_best] = bmc_nm(X_miss_random,[1 1], ...
row_laplacian,col_laplacian,'approximate',true,'sample_size',10,'show',false);
imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum BIC is %d .\n',bic_best);
% AIC as the objective
[X_complete,gamma_best,aic_best] = bmc_nm(X_miss_random,[1 1], ...
row_laplacian,col_laplacian,'approximate',true,'sample_size',10,'show',false,...
'criterion', 'AIC');
%imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum AIC is %d .\n',aic_best);
%% Cross Validation grid search
% BIC as objective
row_gamma_list = exp(linspace(-6, 1, 10));
col_gamma_list = exp(linspace(-6, 1, 10));
[X_complete,gamma_best,bic_best] = bmc_cv(X_miss_random, ...
    row_gamma_list,col_gamma_list,row_laplacian,col_laplacian,'K',5);
%imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum BIC is %d .\n',bic_best);
% AIC as objective
[X_complete,gamma_best,aic_best] = bmc_cv(X_miss_random, ...
    row_gamma_list,col_gamma_list,row_laplacian,col_laplacian,'K',5,...
    'criterion', 'AIC');
imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum AIC is %d .\n',aic_best);
%% Grid search with exact solve
% BIC as objective function
row_gamma_list = exp(linspace(-6, 1, 10));
col_gamma_list = exp(linspace(-6, 1, 10));
[X_complete,gamma_best,bic_best] = bmc_grid(X_miss_random, ...
    row_gamma_list,col_gamma_list,row_laplacian,col_laplacian, ...
    'approximate',false);
%imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum BIC is %d .\n',bic_best);
% AIC as objective function
[X_complete,gamma_best,aic_best] = bmc_grid(X_miss_random, ...
    row_gamma_list,col_gamma_list,row_laplacian,col_laplacian, ...
    'approximate',false,'criterion', 'AIC');
imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum AIC is %d .\n',aic_best);
%% Grid search with Hutchinson estimator
row_gamma_list = exp(linspace(-6, 1, 10));
col_gamma_list = exp(linspace(-6, 1, 10));
% BIC objective
[X_complete,gamma_best,bic_best] = bmc_grid(X_miss_random, ...
    row_gamma_list,col_gamma_list,row_laplacian,col_laplacian, ...
    'approximate',true,'sample_size', 1000);
%imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum BIC is %d .\n',bic_best);
% AIC objective
[X_complete,gamma_best,aic_best] = bmc_grid(X_miss_random, ...
    row_gamma_list,col_gamma_list,row_laplacian,col_laplacian, ...
    'approximate',true,'sample_size', 1000, 'criterion', 'AIC');
imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum AIC is %d .\n',aic_best);
%% Quasi-Newton with exact solve
% Note: initial is [0 0] b/c here is input is log of gamma
% BIC objective 
[X_complete,gamma_best,bic_best] = bmc_qn(X_miss_random,[0 0],...
row_laplacian,col_laplacian,'approximate',false,'show',false);
%imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum BIC is %d .\n',bic_best);
% AIC objective
[X_complete,gamma_best,aic_best] = bmc_qn(X_miss_random,[0 0],...
row_laplacian,col_laplacian,'approximate',false,'show',false,...
'criterion', 'AIC');
imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum AIC is %d .\n',aic_best);
%% Quasi-Newton with Hutchinson approximation
% Note: initial is [0 0] b/c here the input is log of gamma
% BIC objective
[X_complete,gamma_best,bic_best] = bmc_qn(X_miss_random,[0 0],...
row_laplacian,col_laplacian,'approximate',true,'sample_size', 100);
%imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum BIC is %d .\n',bic_best);
% AIC objective
[X_complete,gamma_best,aic_best] = bmc_qn(X_miss_random,[0 0],...
row_laplacian,col_laplacian,'approximate',true,'sample_size', 100, ...
'criterion', 'AIC');
imagesc(X_complete);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum AIC is %d .\n',aic_best);
%% Conjugate Gradient with Hutchinson approximation
% BIC objective
[X_complete,gamma_best,bic_best] = bmc_qn_cg(X_miss_random,[0 0],...
                row_laplacian,col_laplacian,'approximate',true,...
                'sample_size',5,'show', false);
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum BIC is %d .\n',bic_best);
% AIC objective
[X_complete,gamma_best,aic_best] = bmc_qn_cg(X_miss_random,[0 0],...
                row_laplacian,col_laplacian,'approximate',true,...
                'sample_size',5,'show', false, 'criterion', 'AIC');
fprintf('The best gamma will be %d and %d .\n',gamma_best);
fprintf('The mininum BIC is %d .\n',aic_best);