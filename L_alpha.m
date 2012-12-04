function [f, g, h] = L_alpha (x, model, data)

indzero = [];
if(model.phase==1 && model.option==2) %% only for training phase of labeled LDA
    M = [data.annotations];
else
    M = [ones(model.N, model.K)];
end

if(model.phase==1 && model.option==2)
    indzero = find(M==0);  %% for labeled LDA
end

alphaterm  = [repmat(x,model.N,1)];
mod_alpha  = (M.*alphaterm);   %alpha_{ni}  -- we need to ignore co-ordinates which are zero

gmod_alpha = gammaln(mod_alpha);
gmod_alpha(find(gmod_alpha==inf)) = 0;               %gammaln(alpha_{ni}) inf replaced by zero as we need to sum over the co-ordinates
der_alpha  = (M.*alphaterm);

psi_mod_alpha  = psi(mod_alpha);
psi_mod_alpha(find(psi_mod_alpha==-inf)) = 0;

psi_gamma  = psi(model.gamma);
psi_gamma(find(psi_gamma==-inf)) = 0;      % some gammas are zeroed out in training phase, need to undo them here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term1 = sum(gammaln(sum(mod_alpha,2)));
term2 = - sum(sum(gmod_alpha));
term3 = (psi_gamma -repmat(psi(sum(model.gamma,2)),1,model.K)).*(mod_alpha-1);
term3(indzero) = 0;
f     = -(term1+term2+sum(sum(term3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    
    term11 = sum(psi(sum(mod_alpha,2)));
    term12 = -sum(psi_mod_alpha,1);
    term13 = (psi_gamma -repmat(psi(sum(model.gamma,2)),1,model.K));
    term13(indzero) = 0;
    g = -(term11+term12+sum(term13,1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% function hessian %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to do; not yet completed
    if nargout > 2
        
        term21 = sum(psi(1,sum(mod_alpha,2)));
        term22 = -sum(psi(1,mod_alpha),1);
        h = -diag(term21+term22)*(1-model.epsilon)^2;
        
    end
    
end

end
