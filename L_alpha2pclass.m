function [f, g, h] = L_alpha2pclass (x, model, data, classnum)


trind = find(data.classlabels==classnum);
N = length(trind);
gammatrind = model.gamma(trind,:);
x = x(model.k1+1:end);

if(model.phase==1)
    M = [data.annotations(trind,:) ones(N, model.k2)];
else
    M = [ones(N, model.K)];
end

alphaterm1 = [repmat(model.alpha1,N,1) zeros(N, model.k2)];
alphaterm2 = [zeros(N, model.k1) repmat(x,N,1)];
mod_alpha  = (M.*alphaterm1)*model.epsilon + (1-model.epsilon)*alphaterm2;   %alpha_{ni}  -- we need to ignore co-ordinates which are zero

gmod_alpha = gammaln(mod_alpha);
gmod_alpha(find(gmod_alpha==inf)) = 0;               %gammaln(alpha_{ni}) inf replaced by zero as we need to sum over the co-ordinates
der_alpha  = (M.*alphaterm1) - alphaterm2;

psi_mod_alpha  = psi(1,mod_alpha);
psi_mod_alpha(find(psi_mod_alpha==inf)) = 0;

psi_gamma  = psi(1,gammatrind);
psi_gamma(find(psi_gamma==inf)) = 0;      % some gammas are zeroed out in training phase, need to undo them here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term1 = sum(gammaln(sum(mod_alpha,2)));
term2 = - sum(sum(gmod_alpha));
term3 = sum(sum((psi_gamma -repmat(psi(1,sum(gammatrind,2)),1,model.K)).*(mod_alpha-1)));
f     = -(term1+term2+term3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    
    term11 = sum(psi(1,sum(mod_alpha,2)));
    term12 = -sum(psi_mod_alpha,1);
    term13 = [sum(psi_gamma,1) - sum(psi(1,sum(gammatrind,2)))]*(1-model.epsilon);
    g = -(term11+term12+term13);
    g(1:model.k1) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% function hessian %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to do; not yet complete but does not matter
    if nargout > 2
        
        term21 = sum(psi(2,sum(mod_alpha,2)));
        term22 = -sum(psi(2,mod_alpha),1);
        h = -diag(term21+term22)*(1-model.epsilon)^2;
        
    end
    
end

end
