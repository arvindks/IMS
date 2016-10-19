function [X_complete,gamma_best,f_best] = bmc_nm(X_miss,gamma0,...
    row_laplacian,col_laplacian,varargin)
% BMC_NM A function that performs biclustered matrix completion 
%        with tuning parameters selected by minmizing BIC or AIC using 
%        neldermead algorithm
%
% INPUT:
%   X_miss        - n-by-p complete data matrix with missing values indicated 
%                   by NaN
%   gamma0        - vector of initial regularization parameters 
%   row_laplacian - Laplacian matrix for rows
%   col_laplacian - Laplacian matrix for columns
%
% Optional inputs and their default values are:
%   criterion     - Selection criterion "BIC" or "AIC" {"BIC"}
%   approximate   - logical, use Huchinson if true {false} ; 
%   sample_size   - number of samples when using Hutchsinson {1e3};
%   max_iter      - maximum number of iterations when using fminsearch {1e5};
%   tol           - tolerance when using fminsearch {1e-6}   
%   show          - controls the iteration log {false}
%
% Output:
%   X_complete    - n-by-p data matrix with missing entries computed
%   gamma_best    - vector of best gamma
%   f_best        - best objective value
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
params.addParameter('max_iter',400, @(x) x >= 1);
params.addParameter('tol',1e-6, @(x) x > 0);
params.addParameter('show',false);
params.KeepUnmatched = true;
params.parse(varargin{:});

criterion = params.Results.criterion;
approximate = params.Results.approximate;
sample_size = params.Results.sample_size;
max_iter = params.Results.max_iter;
tol = params.Results.tol;
show = params.Results.show;

[n, p] = size(X_miss);
np = n * p;
cx = X_miss(:);
miss_idx = isnan(cx);
numObs = sum(~miss_idx);
P = spdiags(~miss_idx(:),0,np,np);
%P = diag(1 - miss_idx);
Px = cx;
Px(miss_idx) = 0;
Px = sparse(Px);
%LrI = kron(eye(p),row_laplacian);
%LcI = kron(col_laplacian, eye(n));
LrI = kron(speye(p),row_laplacian);
LcI = kron(col_laplacian, speye(n)); 
        
if approximate
    W = 2 * binornd(1, 0.5, np, sample_size) - 1;
end

    function f = cost(gamma)
        gamma1 = gamma(1); gamma2 = gamma(2);
        % capping the gamma below 10^20
        if gamma1 > 1e20
            gamma1 = 1e20;
            warning('Gamma_r has reached 1e20');
        end
        if gamma2 > 1e20
            gamma2 = 1e20;
            warning('Gamma_c has reached 1e20');
        end        
        %Construct the matrix explicitly
        S = P + gamma1 * LrI + gamma2 * LcI;
        
        % use Huchinson to estimate trace 
        if approximate 
            if any(gamma(:) < 1e-16)
                f = 1e20;
            else
                SinvW = S \ W;
                df = mean(sum(W .* SinvW));
                z = S \ Px;
                r = ~miss_idx .* z - Px;   %Compute the residual
                %bic = log(np) * df + sum((z(~miss_idx) - cx(~miss_idx)).^2);
                if strcmp(criterion,'BIC')
                    f = numObs * log(norm(r).^2) + log(numObs) * df; 
                elseif strcmp(criterion,'AIC')
                    f = numObs * log(norm(r).^2) + 2 * df;
                else
                    error('Error. \nOnly BIC or AIC can be selected as Criterion in current implementation ')
                end
            end
        % exact solve    
        else
            if any(gamma(:) < 1e-16)
                f = 1e20;
            else
                Sinv = inv(S);
                df = trace(Sinv);
                z = Sinv * Px;
                r = ~miss_idx .* z - Px;   %Compute the residual
                %bic = log(np) * df  + sum((z(~miss_idx) - cx(~miss_idx)).^2);
                if strcmp(criterion,'BIC')
                    f = numObs * log(norm(r).^2) + log(numObs) * df; 
                elseif strcmp(criterion,'AIC')
                    f = numObs * log(norm(r).^2) + 2 * df;
                else
                    error('Error. \nOnly BIC or AIC can be selected as Criterion in current implementation ')
                end
             end
        end

        if show
            fprintf('\n')
            fprintf('\n Current parameter: Gamma_r =  %f, Gamma_c = %f',...
                    gamma1, gamma2)
            if any(gamma(:) < 1e-16)
                fprintf('\n No df since gamma < 0')
            else
                fprintf('\n Current df = %f', full(df))
            end
            
            fprintf('\n Current objective = %f', full(f))
            fprintf('\n')
        end
    end

[gamma_best, f_best] = fminsearch(@cost, gamma0, ...
    optimset('TolFun',tol,'MaxIter', max_iter));

% calcualte complete matrix using the selected gamma
S = P + gamma_best(1) * kron(eye(p),row_laplacian) + ...
            gamma_best(2) * kron(col_laplacian, eye(n));
z = S \ Px;
X_complete = reshape(z, [n p]);
if show
    fprintf('\n Final parameter: Gamma_r =  %f, Gamma_c = %f',...
                    gamma_best(1), gamma_best(2))
    fprintf('\n Final objective = %f', full(f_best))                
    fprintf('\n')
end
end

