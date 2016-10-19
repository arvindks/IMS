function [lb,ub] = confidence_interval(gamma,x,Po,LcI,LrI,Omega,alpha)
    %% Compute the Confidence interval for the BIC exactly
    % Also computes the function and gradient as a bonus
    % Written by A. Saibaba  
    
    
    gamma1 = gamma(1);  gamma2 = gamma(2);
    
    %Construct the matrix explicitly
    Pox = Po*x;
    S = Po + gamma1*LcI + gamma2*LrI;
    z = S\(Pox);
   
    r = Po*(z-x);   %Compute the residual
    Por = Po*r;
     
    %Contribution from the regularization
    SinvLc = S\(LcI);   SinvLr = S\LrI;
    Sc = SinvLc/S;      Sr = SinvLr/S;
    
    
    %Hessian - trace contribution
    Scc = 2*SinvLc*Sc;
    Srr = 2*SinvLr*Sr;
    Scr = SinvLc*Sr + SinvLr*Sc;
    
    %Contribution to the Hessian from the trace
    Ht  = [trace(Scc) trace(Scr); trace(Scr) trace(Srr)];
    
    %Contributions to the Hessian from the misfit
    ScPox = Sc*Pox; SrPox = Sr*Pox;
    Hm = 2*[Pox'*Scc*Por + ScPox'*Po*ScPox, ...
                Pox'*(Scr*Por) + ScPox'*Po*SrPox; ...
                Pox'*(Scr*Por) + SrPox'*Po*ScPox, ...
                Pox'*(Srr*Por) + SrPox'*Po*SrPox];
            
    %Empirical Hessian
    Jk = zeros(2);      Sigmak = zeros(2);
    np = size(Pox,1);
    lognp = log(np);    ns = size(Omega,2);
    for i = 1:ns
        w = Omega(:,i);
        Jk = Jk +  Hm + lognp*[w'*Scc*w w'*Scr*w; w'*Scr*w w'*Srr*w];
        
        gk = [-4*Pox'*Sc*Por - lognp*(w'*Sc*w); ...
                   -4*Pox'*Sr*Por - lognp*(w'*Sr*w)];
        Sigmak = Sigmak + gk*gk';       
    end
    
    Jk = Jk/ns;
    Sigmak = Sigmak/ns;
    
    V = Jk\(Sigmak/Jk);
    
    %Sum the contribution to the Hessian        
    H = Hm + Ht;   
    
    %Compute eigenvalues
    l   = eig(sqrtm(V)*H*sqrtm(V));
    
    
    %Simulate the confidence intervals
    qs = weighted_chisquared(l,alpha);
    lb = qs(1)/ns;
    ub = qs(2)/ns;
end

function qs= weighted_chisquared(l,alpha)
    ns = 1000;
    samples = l(1)*randn(ns,1).^2 + l(2)*randn(ns,1).^2;
    
    qs = quantile(samples, [alpha/2,1-alpha/2]);
    
end
