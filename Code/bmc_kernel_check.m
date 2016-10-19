function [invertable] = bmc_kernel_check(X_miss,row_laplacian,col_laplacian)
% BMC_KERNEL_CHECK   given X_miss, row and column Laplacian, check if S is 
%                   invertable
%
% NOTE:             Use the "networkComponents.m" written by 
%                   Daniel Larremore. For details, see 
%                   http://www.mathworks.com/matlabcentral/fileexchange/
%                   42040-find-network-components
% INPUT:
%   X_miss        - n-by-p complete data matrix with missing values indicated 
%                   by NaN
%   row_laplacian - Laplacian matrix for rows
%   col_laplacian - Laplacian matrix for columns
%
% Output:
%   invertable    - logical value indicates whether S is invertable or
%                    not
%
% AUTHOR: Liuyi Hu (lhu@ncsu.edu)
%
cx = X_miss(:);
miss_idx = isnan(cx);
p = ~miss_idx;
% Convert the Laplacian matrix to adjacency matrix
A_r = (row_laplacian~=0);
N_r = size(A_r,1);
A_r(1:N_r+1:end) = 0;
A_c = (col_laplacian~=0);
N_c = size(A_c,1);
A_c(1:N_c+1:end) = 0;

% get the connected components of rows and columns
[nComponents_r,sizes_r,members_r] = networkComponents(A_r);
[nComponents_c,sizes_c,members_c] = networkComponents(A_c);
% test if P_{Omega}[X_{Bc} \kron X_{Ar}] = 0
invertable = true;
for row_idx = 1:nComponents_r
    for col_idx = 1:nComponents_c
        row_ind_vec = zeros(1, N_r);
        row_ind_vec(members_r{row_idx})  = 1;
        col_ind_vec = zeros(1, N_c);
        col_ind_vec(members_c{col_idx})  = 1;   
        Ker_rc = kron(col_ind_vec,row_ind_vec)';
        if norm(Ker_rc .* p) == 0
            invertable = false;
            break
        end
    end
    if ~invertable
        break
    end
end

if ~invertable
  warning('S is not invertible')
  fprintf('\n');
else
  fprintf('S is invertible \n');
end

end
