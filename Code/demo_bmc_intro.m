%% This is an real data analysis
%  Data: This is the Mondrian data
%        which is a 370-by-380 matrix.
%  This file compares matrix completion via singular value thresholding and
%  biclustered matrix completion
%              
%           Comparison metrics inludes:
%           a. MSE over missing entries
%
%  Date: 2016-05-03
%  Author: Eric Chi eric_chi@ncsu.edu
%

clear;
close all;
rho_tol = 1e-4;

filename = '/Users/ecchi/Dropbox/Work/Research/00_Active/00_NCSU/01_BMC/Data/';
A = imread([filename,'compositionabypietmondrian.jpg']);
B = rgb2gray(A);
X = double(B);
X_true = X;

filename = '/Users/ecchi/Dropbox/Work/Research/00_Active/00_NCSU/01_BMC/Data/Mondrian/';
% Laplacian matrix
row_laplacian = dlmread([filename,'Lr_5.txt'],',', 1, 1);
col_laplacian = dlmread([filename,'Lc_5.txt'],',', 1, 1);

n_repeat = 50;
rho = 0.5;
gridpts = 11;

mse_svt = zeros(n_repeat,1);
mse_bmc = zeros(n_repeat,1);
sigma = 0.5*255;
for repeat_idx = 1:n_repeat
    repeat_idx
    tic;
    %% Add noise
    X = X_true + sigma*randn(size(X_true));
    
    %% Mask missing entries
    tol = 1e-4;
    [n, p] = size(X);
    missfrac = 0.5;
    valfrac = 0.1;
    imiss = randsample(n*p, round(n*p*missfrac), false);
    nMissing = length(imiss);
    ival   = randsample(setdiff(1:n*p, imiss), round(n*p*valfrac), false);

    X_validate = X(ival);
    X_miss = X;
    X_miss(union(ival,imiss)) = NaN;
    
    %% find max lambda for nuclear norm regularization
    [dummy,stats] = matrix_impute_Nesterov(X_miss,inf);
    maxlambda = stats.maxlambda;

%% optimization on grid pts with warm start

    lambdas = maxlambda * (rho.^(0:(gridpts-1)));

    nLambda = length(lambdas);
    miss_rate = zeros(nLambda,1);
    Zi = cell(gridpts,1);
    ranks = zeros(gridpts,1);
    objval = zeros(gridpts,1);
    iterations = zeros(gridpts,1);



    validation_error = zeros(gridpts,1);

    for i=1:gridpts
        lambda = lambdas(i);
        if (i==1)
            Y0 = [];
        else
            Y0 = Z;
        end
        [Z,stats] = matrix_impute_Nesterov(X_miss,lambda,'Y0',Y0,'Display','off','TolFun',tol);
        ranks(i) = stats.rank;
        objval(i) = stats.objval;
        iterations(i) = stats.iterations;
        validation_error(i) = sum( (Z(ival) - X(ival)).^2 );
        Zi{i} = Z;
    end

    idx_min = find(validation_error == min(validation_error));

    % MSE SVT
    mse_svt(repeat_idx) = sum( ( Zi{idx_min}(imiss) - X_true(imiss)).^2 )/length(imiss);
    
    X_miss = X;
    X_miss(imiss) = NaN;
    sample_size = 5;
            [X_complete,gamma_best,bic_best] = bmc_qn(X_miss,[0 0],...
                row_laplacian,col_laplacian,'approximate',true,...
                'sample_size',sample_size,'show', false);
            
    % MSE BMC
    mse_bmc(repeat_idx) = sum( ( X_complete(imiss) - X_true(imiss)).^2 )/length(imiss);
    toc
end

csvwrite('mondrian_miss.csv',X_miss);
csvwrite('mondrian_svt.csv',Zi{idx_min});

csvwrite('mondrian_mse_svt.csv',mse_svt);
csvwrite('mondrian_mse_bmc.csv',mse_bmc);

X_bmc = zeros(n,p);
for j = 1:n*p
    X_bmc(j) = X_complete(j);
end
csvwrite('mondrian_bmc.csv',X_bmc);

figure;
imagesc(X_true)
colormap gray;
set(gca,'XTick',[]) % Remove the ticks in the x axis!
set(gca,'YTick',[]) % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]) % Make the axes occupy the hole figure
saveas(gcf,'Mondrian_True','pdf')
print -depsc -tiff -r300 Mondrian_True;

filename = '/Users/ecchi/Dropbox/Work/Research/00_Active/00_NCSU/01_BMC/Manuscript/Figures/';
print([filename,'Mondrian_True'],'-depsc2','-r300');

figure;
imagesc(X_miss);
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')
axis off

filename = '/Users/ecchi/Dropbox/Work/Research/00_Active/00_NCSU/01_BMC/Manuscript/Figures/';
print([filename,'Mondrian_Missing_Noisy'],'-depsc2','-r300');


figure;
imagesc(Zi{idx_min})
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')
axis off
 
filename = '/Users/ecchi/Dropbox/Work/Research/00_Active/00_NCSU/01_BMC/Manuscript/Figures/';
print([filename,'Mondrian_SVT'],'-depsc2','-r300');


figure;
imagesc(X_complete)
axis
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')
axis off

filename = '/Users/ecchi/Dropbox/Work/Research/00_Active/00_NCSU/01_BMC/Manuscript/Figures/';
print([filename,'Mondrian_BMC'],'-depsc2','-r300');


plot(X_true(imiss),Zi{idx_min}(imiss),'.r')
plot(X_true(imiss), X_complete(imiss),'.r')

% Visualize prediction error
plot( Zi{idx_min}(imiss) - X(imiss), '.r')
figure; plot( X_complete(imiss) - X(imiss), '.r')

plot(X_complete(imiss), Zi{idx_min}(imiss), '.r')
