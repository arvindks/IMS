function [f,grad,H] = hessian_exact(gamma,x,Po,LcI,LrI)
    %% Compute the Hessian of the BIC exactly
    % Also computes the function and gradient as a bonus
    % Written by A. Saibaba  
    gamma1 = gamma(1);  gamma2 = gamma(2);
    
    %Number of samples
    ns = size(Omega,2);
    
    %Size of problem
    np = size(Po,1);
    
    %Construct the matrix explicitly
    Pox = Po*x;
    S = Po + gamma1*LcI + gamma2*LrI;
    z = S\(Pox);
    
    
    r = Po*(z-x);   %Compute the residual
    f = norm(r).^2 + trace(inv(S)); %Compute objective function

    %Contribution from the misfit
    Por = Po*r;
    SinvPor = S\(Por); 
    drdg1 = -2*z'*(LcI*SinvPor);
    drdg2 = -2*z'*(LrI*SinvPor);
    
    %Contribution from the regularization
    SinvLc = S\(LcI);   SinvLr = S\LrI;
    Sc = SinvLc/S;      Sr = SinvLr/S;
    
    dtdg1 = -trace(Sc);
    dtdg2 = -trace(Sr);
    
    %Compute the gradient
    grad = [drdg1 + dtdg1; drdg2 + dtdg2];
    
    Scc = 2*SinvLc*Sc;
    Srr = 2*SinvLr*Sr;
    Scr = SinvLc*Sr + SinvLr*Sc;
    
    %Contribution to the Hessian from the trace
    Ht  = [trace(Scc) trace(Scr); trace(Scr) trace(Srr)];
    
    %Contribution to the Hessian from the misfit
    ScPox = Sc*Pox; SrPox = Sr*Pox;
    Hm = 2*[Pox'*Scc*Por + ScPox'*Po*ScPox, ...
                Pox'*(Scr*Por) + ScPox'*Po*SrPox; ...
                Pox'*(Scr*Por) + SrPox'*Po*ScPox, ...
                Pox'*(Srr*Por) + SrPox'*Po*SrPox];
            
            
    %Sum the contribution to the Hessian        
    H = Hm + Ht;        
            
end