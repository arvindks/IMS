function L = laplacian(n,cluster,varargin)
% LAPLACIAN A function that performs biclustered matrix completion 
%        with tuning parameters selected by minmizing BIC using 
%        neldermead algorithm
%
% INPUT:
%   n - dimension of laplacian matirx
%   cluster - a vector inindicates which columns are within same cluster
%             the number in the cluster should sum up to the dimension 
%             of the  dimension e.g. c(50, 50) means the first 50 columns 
%             are within the same  cluster and the next 50 columns are 
%             within another cluster, in this  case, n should be 50+50=100
%
% Optional inputs and their default values are:
%   weight_in - weight if within the same group, {1}
%   weight_out - weight if not in the same group, {0.001}
%
% Output:
%   L - A n-by-n laplacian matrix
%
% AUTHOR: Liuyi Hu (lhu@ncsu.edu)
params = inputParser;
params.addParameter('weight_in',1, @(x) x > 0);
params.addParameter('weight_out',0.001, @(x) x >= 0);

params.KeepUnmatched = true;
params.parse(varargin{:});

weight_in = params.Results.weight_in;
weight_out = params.Results.weight_out;
% check condition
if sum(cluster) ~= n
    error('Error. \n Sum of cluster vector should equal to dimension')
end

% L = D - W
W = weight_out * ones(n, n);
start_idx = 1;
for i=1:length(cluster)
    end_idx = sum(cluster(1:i));
    W(start_idx:end_idx,start_idx:end_idx) = weight_in;
    start_idx = end_idx + 1;
end
W(1:n+1:end) = 0;
d = sum(W, 1);
L = diag(d) - W;
end

