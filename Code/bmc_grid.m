function [X_complete,gamma_best,f_best,f_grid] = bmc_grid(X_miss,...
    row_gamma_list,col_gamma_list,row_laplacian,col_laplacian,varargin)
% BMC_GRID A function that performs biclustered matrix completion 
%        with tuning parameters selected by minmizing BIC using 
%        grid search
%
% INPUT:
%   X_miss         - n-by-p complete data matrix with missing values  
%                    indicated by NaN
%   row_gamma_list - vector of row clustering regularization parameter 
%   col_gamma_list - vector of column clustering regularization parameter 
%   row_laplacian  - Laplacian matrix for rows
%   col_laplacian  - Laplacian matrix for columns
%
% Optional inputs and their default values are:
%   criterion     - Selection criterion "BIC" or "AIC" {"BIC"}
%   approximate    - logical, use Huchinson if true {false} ; 
%   sample_size    - number of samples when using Hutchsinson {1e3};
%  
% Output:
%   X_complete     - n-by-p data matrix with missing entries computed
%   gamma_best     - vector of best gamma
%   f_best         - best objective value
%   f_grid         - grid of objective values 
%
% Modifications:
%  08/24/2016     - Change the BIC objective function from 
%                   norm(x-z)^2 +log(np)* df
%                   to 
%                   |\Omega|log(norm(x-z)^2) + log(|\Omega|) * df
%                 - Add an option to calcuate AIC
%                                                       
%
% AUTHOR: Liuyi Hu (lhu@ncsu.edu)
%
params = inputParser;
params.addParameter('criterion', 'BIC')
params.addParameter('approximate',false);
params.addParameter('sample_size',1e3, @(x) x >= 1);

params.KeepUnmatched = true;
params.parse(varargin{:});

criterion = params.Results.criterion;
approximate = params.Results.approximate;
sample_size = params.Results.sample_size;

[n, p] = size(X_miss);
np = n * p;
cx = X_miss(:);
miss_idx = isnan(cx);
numObs = sum(~miss_idx);
P = diag(1 - miss_idx);
Px = cx;
Px(miss_idx) = 0;
% initialize to save 
f_grid = zeros(length(row_gamma_list),length(col_gamma_list));

if approximate
    W = 2 * binornd(1, 0.5, np, sample_size) - 1;
end

% loop over row_gamma_list and col_gamma_list
for i = 1:length(row_gamma_list)
    for j=1:length(col_gamma_list)
        
        % construct S matrix
        S = P + row_gamma_list(i) * kron(eye(p),row_laplacian) + ...
            col_gamma_list(j) * kron(col_laplacian, eye(n));
        if approximate
            df =  mean(sum(W .* (S \ W)));
            z = S \ Px;
            r = ~miss_idx .* z - Px;   %Compute the residual
            %bic = log(np) * df + sum((z(~miss_idx) - cx(~miss_idx)).^2);
        else
            Sinv = inv(S);
            df = trace(Sinv);
            z = Sinv * Px;
            r = ~miss_idx .* z - Px;   %Compute the residual
            %bic = log(np) * df  + sum((z(~miss_idx) - cx(~miss_idx)).^2);
        end
        if strcmp(criterion,'BIC')
            f = numObs * log(norm(r).^2) + log(numObs) * df; 
        elseif strcmp(criterion,'AIC')
            f = numObs * log(norm(r).^2) + 2 * df;
        else
            error('Error. \nOnly BIC or AIC can be selected as Criterion in current implementation ')
        end
        f_grid(i, j) = f;
        
        
    end
end
  
    
[f_best,I] = min(f_grid(:));
[I_row, I_col] = ind2sub(size(f_grid),I);
gamma_best = [row_gamma_list(I_row) col_gamma_list(I_col)];

% calcualte complete matrix using the selected gamma
S = P + gamma_best(1) * kron(eye(p),row_laplacian) + ...
            gamma_best(2) * kron(col_laplacian, eye(n));
z = S \ Px;
X_complete = reshape(z, [n p]);

end

