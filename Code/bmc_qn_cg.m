function [X_complete,gamma_best,f_best,output] = bmc_qn_cg(X_miss,loggamma0,...
    row_laplacian,col_laplacian,varargin)
% BMC_QN_CG Solve bmc problem using quasi-newton method with log transform'
%           and use incomplete cholesky and conjuate gradient to deal with
%           large and sparse matrix
% NOTE: the parameter here is the log of gamma!! 
% INPUT:
%   X_miss        - n-by-p complete data matrix with missing values indicated 
%                   by NaN
%   loggamma0     - vector of initial regularization parameters 
%   row_laplacian - Laplacian matrix for rows
%   col_laplacian - Laplacian matrix for columns
%
% Optional inputs and their default values are:
%   criterion     - Selection criterion "BIC" or "AIC" {"BIC"}
%   approximate   - logical, use Hutchinson if true {true} ; 
%   sample_size   - number of samples when using Hutchsinson {1e3};
%   max_iter      - maximum number of iterations when using fminsearch {1e5};
%   tol           - tolerance when using cgLanczos {1e-4}
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
params.addParameter('approximate',true);
params.addParameter('sample_size',1e3, @(x) x >= 1);
params.addParameter('max_iter',1e3, @(x) x >= 1);
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
LrI = kron(speye(p),row_laplacian);
LcI = kron(col_laplacian, speye(n));

if approximate
    W = 2 * binornd(1, 0.5, np, sample_size) - 1;
    max_iter_cg = 10 * np;
end

    function [f,grad] = gradient_exact(loggamma)
    % GRADIENT_EXACT gradient of BMC problem
    % AUTHOR: Arvind Krishna Saibaba asaibab@ncsu.edu
    %         Liuyi Hu lhu@ncsu.edu
        gamma1 = exp(loggamma(1));  gamma2 = exp(loggamma(2)); 
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
        opts.droptol = 1.e-2;
        L = ichol(S,opts);
        
        if strcmp(criterion,'BIC')
            if approximate
            
                fprintf('\n');
                z = pcg(S,Px,tol,max_iter_cg,L,L');
            
                %[z,~] = cgLanczos(S,Px,0,0,max_iter_cg,tol);
                %z = S \ Px;
                r = ~miss_idx .* z - Px; %Compute the residual
            
                SinvW = NaN(np, sample_size);
            
                for sample_idx = 1:sample_size
                    w = W(:,sample_idx);
                    Sinvw = pcg(S,w,tol,max_iter_cg,L,L');
                    SinvW(:,sample_idx) = Sinvw;
                end
 
                df = mean(sum(W .* SinvW));
                %f = norm(r).^2 + log(np) * df;
            
                if nargout > 1 % gradient required
                    %Contribution from the misfit
                    SinvPor =  pcg(S,r,tol,max_iter_cg,L,L');
                    drdg1 = - gamma1*2*(numObs/norm(r).^2)*z'*(LrI*SinvPor);
                    drdg2 = - gamma2*2*(numObs/norm(r).^2)*z'*(LcI*SinvPor);
    
                    %Contribution from the regularization
                    dtdg1 = - gamma1*log(numObs) * mean(sum(SinvW .* (LrI*SinvW)));
                    dtdg2 = - gamma2*log(numObs) * mean(sum(SinvW .* (LcI*SinvW)));
    
                    %Compute the gradient
                    grad = [drdg1 + dtdg1; drdg2 + dtdg2];
                end
            else
                fprintf('\n')
                fprintf('Have only implemented Conjugate Gradient with Hutchinson \n')
            end
            f = numObs * log(norm(r).^2) + log(numObs) * df;
        elseif strcmp(criterion,'AIC')
            if approximate
            
                fprintf('\n');
                z = pcg(S,Px,tol,max_iter_cg,L,L');
            
                %[z,~] = cgLanczos(S,Px,0,0,max_iter_cg,tol);
                %z = S \ Px;
                r = ~miss_idx .* z - Px; %Compute the residual
            
                SinvW = NaN(np, sample_size);
            
                for sample_idx = 1:sample_size
                    w = W(:,sample_idx);
                    Sinvw = pcg(S,w,tol,max_iter_cg,L,L');
                    SinvW(:,sample_idx) = Sinvw;
                end
 
                df = mean(sum(W .* SinvW));
                %f = norm(r).^2 + log(np) * df;
            
                if nargout > 1 % gradient required
                    %Contribution from the misfit
                    SinvPor =  pcg(S,r,tol,max_iter_cg,L,L');
                    drdg1 = - gamma1*2*(numObs/norm(r).^2)*z'*(LrI*SinvPor);
                    drdg2 = - gamma2*2*(numObs/norm(r).^2)*z'*(LcI*SinvPor);
    
                    %Contribution from the regularization
                    dtdg1 = - gamma1*2 * mean(sum(SinvW .* (LrI*SinvW)));
                    dtdg2 = - gamma2*2 * mean(sum(SinvW .* (LcI*SinvW)));
    
                    %Compute the gradient
                    grad = [drdg1 + dtdg1; drdg2 + dtdg2];
                end
            else
                fprintf('\n')
                fprintf('Have only implemented Conjugate Gradient with Hutchinson \n')
            end
            f = numObs * log(norm(r).^2) + 2 * df;
        else
            error('Error. \nOnly BIC or AIC can be selected as Criterion in current implementation ')
        end
        if show
            fprintf('\n')
            fprintf('\n Current parameter: Gamma_r =  %f, Gamma_c = %f',...
                    gamma1, gamma2)
            fprintf('\n Current df = %f', full(df))
            fprintf('\n Current objective = %f', full(f))
            fprintf('\n')
        end
        
    end

options = optimoptions('fminunc', 'Algorithm','quasi-newton',...
    'Gradobj','on','MaxIter',max_iter,'TolFun',tol);
[loggamma_best, f_best,~,output] = fminunc(@gradient_exact, loggamma0, options);

gamma_best = exp(loggamma_best);
% calcualte complete matrix using the selected gamma
S = P + gamma_best(1) * LrI + gamma_best(2) * LcI;
%[z,~] = cgLanczos(S,Px,0,0,max_iter_cg,tol);
opts.droptol = 1.e-2;
L = ichol(S,opts);
z = pcg(S,Px,tol,max_iter_cg,L,L');
X_complete = reshape(z, [n p]);

if show
    fprintf('\n Final parameter: Gamma_r =  %f, Gamma_c = %f',...
                    gamma_best(1), gamma_best(2))
    fprintf('\n Final objective = %f', full(f_best))                
    fprintf('\n')
end

end


