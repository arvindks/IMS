function [X_complete,gamma_best,f_best,mse_grid] = bmc_cv(X_miss,...
    row_gamma_list,col_gamma_list,row_laplacian,col_laplacian,varargin)
% BMC_CV A function that performs biclustered matrix completion 
%        with tuning parameters selected by K-fold cross validation
%        
%
% INPUT:
%   X_miss         - n-by-p complete data matrix with missing values indicated 
%                    by NaN
%   row_gamma_list - vector of row clustering regularization parameter 
%   col_gamma_list - vector of column clustering regularization parameter 
%   row_laplacian  - Laplacian matrix for rows
%   col_laplacian  - Laplacian matrix for columns
%
% Optional inputs and their default values are:
%   criterion     - Selection criterion "BIC" or "AIC" {"BIC"}
%   K              - the number of folds in cross validation {5}
%   seed           - seed for the K-fold partition
%  
% Output:
%   X_complete     - n-by-p data matrix with missing entries computed
%   gamma_best     - vector of best gamma
%   f_best         - best objective value
%   mse_grid       - grid of MSE values 
%
% Modifications:
%  08/24/2016     - Change the BIC objective function from 
%                   norm(x-z)^2 +log(np)* df
%                   to 
%                   |\Omega|log(norm(x-z)^2) + log(|\Omega|) * df
%                 - Add an option to calcuate AIC
%                                                       
%
% Date: 2016-04-23
% AUTHOR: Liuyi Hu (lhu@ncsu.edu)
%
params = inputParser;
params.addParameter('criterion', 'BIC')
params.addParameter('K',5, @(x) x >= 1)
params.addParameter('seed',1991)
params.KeepUnmatched = true;
params.parse(varargin{:});

criterion = params.Results.criterion;
K = params.Results.K;
seed = params.Results.seed;

[n, p] = size(X_miss);
np = n * p;
cx = X_miss(:);
miss_idx = isnan(cx);
numObs = sum(~miss_idx);
% K fold cv partition
obs_Indices = find(miss_idx ~= 1);
rng(seed)
CVO = cvpartition(length(obs_Indices),'k',K);
cv_mse = zeros(CVO.NumTestSets,1);

% initialize to save 
mse_grid = zeros(length(row_gamma_list),length(col_gamma_list));

% loop over row_gamma_list and col_gamma_list
for i = 1:length(row_gamma_list)
    for j=1:length(col_gamma_list)
        
        S0 = row_gamma_list(i) * kron(eye(p),row_laplacian) + ...
             col_gamma_list(j) * kron(col_laplacian, eye(n));
        
        for l = 1:CVO.NumTestSets
            
            teIdx = CVO.test(l);
            cv_obs_idx = ~miss_idx;
            cv_obs_idx(obs_Indices(teIdx))= 0;
            P = diag(cv_obs_idx);
            Px = cx;
            Px(miss_idx) = 0;
            % set test data as missing
            Px(obs_Indices(teIdx))= 0;
            % construct S matrix
            S = P + S0;
            z = S \ Px;
            % evaluate MSE on test data
            cv_mse(l) = norm(z(obs_Indices(teIdx)) - ...
                        cx(obs_Indices(teIdx))).^2/sum(miss_idx);
            
            
        end
        
        mse_grid(i, j) = mean(cv_mse);

    end
end
  
 fprintf('min mse grid %f \n',min(mse_grid(:)));   
[~,I] = min(mse_grid(:));
[I_row, I_col] = ind2sub(size(mse_grid),I);
gamma_best = [row_gamma_list(I_row) col_gamma_list(I_col)];

% calcualte complete matrix using the selected gamma
P = diag(1 - miss_idx);
S = P + gamma_best(1) * kron(eye(p),row_laplacian) + ...
            gamma_best(2) * kron(col_laplacian, eye(n));
Sinv = inv(S);
df = trace(Sinv);
z = Sinv * Px;
r = ~miss_idx .* z - Px;   %Compute the residual
if strcmp(criterion,'BIC')
    f_best = numObs * log(norm(r).^2) + log(numObs) * df; 
elseif strcmp(criterion,'AIC')
    f_best = numObs * log(norm(r).^2) + 2 * df;
else
    error('Error. \nOnly BIC or AIC can be selected as Criterion in current implementation ')
end
X_complete = reshape(z, [n p]);

end

