function X_miss = missing_pattern(X_complete,missfrac,varargin)
% MISSINGPATTERN Convert a complete matrix to a matrix with missing 
% cells following some pattern  
%
% INPUT:
%   X_complete - n-by-p complete data matrice
%   missfrac - fraction of missing entries
%   
% Optional inputs and their default values are:
%   pattern -  "random" - missing at random; 
%              "row" - missing entire row at random; 
%              "both" - missing entire row as well as missing at random.
%              {"random"}
%   
%             
% Output:
%   X_miss - n-by-p data matrice with missing entries denoted as NaN
%
% AUTHOR: Liuyi Hu (lhu@ncsu.edu)

params = inputParser;
params.addParameter('pattern','random');
params.KeepUnmatched = true;
params.parse(varargin{:});

pattern = params.Results.pattern;
[n, p] = size(X_complete);
if strcmp(pattern,'random')
    imiss = randsample(n*p, round(n*p*missfrac), false);
    X_miss = X_complete();
    X_miss(imiss) = NaN;
    X_miss = reshape(X_miss,n,p);
    
elseif strcmp(pattern,'row')
    imissrow = randsample(n, round(n*missfrac), false);
    X_miss = X_complete;
    X_miss(imissrow,:) = NaN;
    
else
    % if none of the conditions is true '
    fprintf('"pattern" should be "random" or "row" \n');

end

