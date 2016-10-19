function [f,grad] = gradient_exact(gamma,x,Po,LcI,LrI)
% NOTE: This function is no longer useful
% 
    gamma1 = gamma(1);  gamma2 = gamma(2);
    miss_idx = isnan(x);
    numObs = sum(~miss_idx);
    %Construct the matrix explicitly
    S = Po + gamma1*LcI + gamma2*LrI;
    z = S\(Po*x);
    
    
    r = Po*(z-x);   %Compute the residual
    f = norm(r).^2 + trace(inv(S)); %Compute objective function

    %Contribution from the misfit
    SinvPor = S\(Po*r); 
    drdg1 = -2*numObs*z'*(LcI*SinvPor)/norm(r).^2;
    drdg2 = -2*numObs*z'*(LrI*SinvPor)/norm(r).^2;
    
    %Contribution from the regularization
    dtdg1 = -log(numObs) * trace(S\(LcI/S));
    dtdg2 = -log(numObs) * trace(S\(LrI/S));
    
    %Compute the gradient
    grad = [drdg1 + dtdg1; drdg2 + dtdg2];
end