function [history,searchdir] = fminunc_path(X_miss,loggamma0,...
    row_laplacian,col_laplacian,varargin)
% FMINUNC_PATH Save the values at each iteration when using fminunc to plot
%               the search path
% INPUT:
%   X_miss - n-by-p complete data matrix with missing values indicated 
%            by NaN
%   loggamma0 - vector of initial regularization parameters 
%   row_laplacian - Laplacian matrix for rows
%   col_laplacian - Laplacian matrix for columns
%
% Optional inputs and their default values are:
%   approximate - logical, use Huchinson if true {false} ; 
%   sample_size - number of samples when using Hutchsinson {1e3};
%   max_iter - maximum number of iterations when using fminsearch {1e5};
%   tol - tolerance when using cgLanczos {1e-4}         
%
% Output:
%   history - The output history is a structure that contains two fields:
%             'x' and 'fval'. 
%             The 'x' field of history contains the sequence of points
%             computed by the algorithm. 
%             The 'fval' field contains the objective function values
%             corresponding to the sequence of points computed.
%   searchdir - the search directions for fminunc at each iteration
%
% AUTHOR: Liuyi Hu (lhu@ncsu.edu)
% Change: 
%         1. Move W out of "gradient_exact" so that at each iteration
%            of quasi-newton, they use the same W (03/21/2016)
%
%

params = inputParser;
params.addParameter('approximate',false);
params.addParameter('sample_size',1e3, @(x) x >= 1);
params.addParameter('max_iter',1e3, @(x) x >= 1);
params.addParameter('tol',1e-6, @(x) x > 0);
params.KeepUnmatched = true;
params.parse(varargin{:});

approximate = params.Results.approximate;
sample_size = params.Results.sample_size;
max_iter = params.Results.max_iter;
tol = params.Results.tol;

[n, p] = size(X_miss);
np = n * p;
cx = X_miss(:);
miss_idx = isnan(cx);
P = diag(1 - miss_idx);
Px = cx;
Px(miss_idx) = 0;
W = 2 * binornd(1, 0.5, np, sample_size) - 1;

    function [f,grad] = gradient_exact(loggamma)
    % GRADIENT_EXACT gradient of BMC problem
    % AUTHOR: Arvind Krishna Saibaba asaibab@ncsu.edu
    %         Liuyi Hu lhu@ncsu.edu
        gamma1 = exp(loggamma(1));  gamma2 = exp(loggamma(2));
        LrI = kron(eye(p),row_laplacian);
        LcI = kron(col_laplacian, eye(n));
        
        %Construct the matrix explicitly
        S = P + gamma1 * LrI + gamma2 * LcI;
        
        if approximate
 
            z = S \ Px;
            r = P * z - Px;   %Compute the residual
            %W = 2 * binornd(1, 0.5, np, sample_size) - 1;
            Pr = P * r;
            Augmented = horzcat(Pr,W);
            SinvAug = S \ Augmented;
            SinvPor = SinvAug(:,1);
            SinvW = SinvAug(:,2:end);
            df = mean(sum(W .* SinvW));
            f = norm(r).^2 + log(np) * df;
            
            if nargout > 1 % gradient required
                %Contribution from the misfit
          
                drdg1 = - gamma1*2*z'*(LrI*SinvPor);
                drdg2 = - gamma2*2*z'*(LcI*SinvPor);
    
                %Contribution from the regularization
                dtdg1 = - gamma1*log(np) * mean(sum(SinvW .* (LrI*SinvW)));
                dtdg2 = - gamma2*log(np) * mean(sum(SinvW .* (LcI*SinvW)));
    
                %Compute the gradient
                grad = [drdg1 + dtdg1; drdg2 + dtdg2];
            end
        else
            z = S \ Px;
            r = P * z - Px;   %Compute the residual
            f = norm(r).^2 + log(np) * trace(inv(S)); %Compute objective function
        
            if nargout > 1 % gradient required
                %Contribution from the misfit
                SinvPor = S \ (P * r); 
                drdg1 = - gamma1*2*z'*(LrI*SinvPor);
                drdg2 = - gamma2*2*z'*(LcI*SinvPor);
    
                %Contribution from the regularization
                dtdg1 = - gamma1*log(np) * trace(S\(LrI/S));
                dtdg2 = - gamma2*log(np) * trace(S\(LcI/S));
    
                %Compute the gradient
                grad = [drdg1 + dtdg1; drdg2 + dtdg2];
            end
            
        end        
        
    end
% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
searchdir = [];
% call optimization
options = optimoptions('fminunc','OutputFcn',@outfun, 'Algorithm',...
    'quasi-newton','Gradobj','on','MaxIter',max_iter,'TolFun',tol);

[~, ~] = fminunc(@gradient_exact, loggamma0, options);
 
 function stop = outfun(x,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
         % Concatenate current search direction with 
         % searchdir.
           searchdir = [searchdir;... 
                        optimValues.searchdirection'];
           plot(x(1),x(2),'o');
         % Label points with iteration number and add title.
         % Add .15 to x(1) to separate label from plotted 'o'
           text(x(1)+.15,x(2),... 
                num2str(optimValues.iteration));
           title('Sequence of Points Computed by fminunc');
         case 'done'
             hold off
         otherwise
     end
 end
 

end


