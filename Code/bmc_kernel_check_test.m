%% This is to test if "bmc_kernel_check.m" work
%  Date: 2016-05-15
%  Author: Liuyi Hu lhu@ncsu.edu
%
%% Set up the missing matrix 
blocks = [25 25 25 25];
X11 = -25 * ones(blocks(1), blocks(1)); %1st block
X12 = 10 * ones(blocks(2), blocks(2));  %2nd blcok
X21 = -10 * ones(blocks(3), blocks(3)); %3rd block
X22 = 25 * ones(blocks(4), blocks(4));  %4th block
X = [X11 X12; X21 X22];
[n, p] = size(X);
% add some noise
X_noise = X + normrnd(0,1,n,p);
% Laplacian
row_laplacian = laplacian(n, [blocks(1), blocks(2)],'weight_out', 0);
col_laplacian = laplacian(p, [blocks(3), blocks(4)],'weight_out', 0);
%% if bicluster induced by the row and column graph Laplacians is 
% completely missing
X_miss_1 = X_noise;
X_miss_1(1:25,1:25) = NaN;
[invertable] = bmc_kernel_check(X_miss_1,row_laplacian,col_laplacian);

%% if bicluster induced by the row and column graph Laplacians is NOT 
% completely missing
X_miss_2 = X_noise;
X_miss_2(1:25,1:25) = NaN;
X_miss_2(25,25) = 1.0;
[invertable] = bmc_kernel_check(X_miss_2,row_laplacian,col_laplacian);

