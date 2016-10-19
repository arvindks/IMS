clc
close all

%Test the Hessian and confidence interval computation


%% Set up the missing matrix 
blocks = [25 25 25 25];
X11 = -25 * ones(blocks(1), blocks(1)); %1st block
X12 = 10 * ones(blocks(2), blocks(2));  %2nd blcok
X21 = -10 * ones(blocks(3), blocks(3)); %3rd block
X22 = 25 * ones(blocks(4), blocks(4));  %4th block
X = [X11 X12; X21 X22];

[n,p] = size(X);

% add some noise
X_noise = X + normrnd(0,1,n,p);

% Laplacian
row_laplacian = laplacian(n, [blocks(1), blocks(2)]);
col_laplacian = laplacian(p, [blocks(3), blocks(4)]);

LrI = kron(speye(p),row_laplacian);
LcI = kron(col_laplacian, speye(n));

%Generate missingness pattern
X_miss= missing_pattern(X_noise, 0.2, 'pattern','random');
x = X_noise(:);

cx = X_miss(:);     miss_idx = isnan(cx);
Px = cx;    Px(miss_idx) = 0;   Px = spdiags([Px],1,n*p,n*p);


%%Optimize BIC to obtain gamma's
% Inputs to the optimiztion problem
gamma0 = [1 1];
loggamma0 = [0 0];
missfrac = 0.2;
n_repeat = 10;
tic;
[X_complete,gamma_qn,bic_best] = bmc_qn(X_miss,loggamma0,...
                            row_laplacian,col_laplacian, 'approximate',true,'sample_size', 5);
z = X_complete(:);



%% OBtain confidence interval for the BIC
ns = 5; gamma = gamma_qn;   alpha = 0.05;   % 95\% confidene interval 
Omega = randn(n*p,ns);
[lb,ub] = confidence_interval(gamma,x,Px,LcI,LrI,Omega,alpha);
