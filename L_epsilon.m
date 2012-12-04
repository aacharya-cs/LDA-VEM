function [f, g, h] = L_epsilon (x, model, data)

indzero = [];

if(model.phase==1)
    M = [data.annotations ones(model.N, model.K-model.k1)];
else
    M = [ones(model.N, model.K)];
end

if(model.phase==1 && model.option>=4)
    indzero = find(M==0);
end

alphaterm1 = [repmat(model.alpha1,model.N,1) zeros(model.N, model.K-model.k1)];
alphaterm2 = [zeros(model.N, model.k1) repmat(model.alpha2,model.N,1)];

if(model.phase==1 && model.option==5)
    ind2  = ones(model.N,model.k1);
    ind2  = [ind2 2*ones(model.N,model.K-model.k1)];
    ind21 = repmat(data.classlabels-1,1,(model.k2/data.Y))*(model.k2/data.Y);
    ind22 = repmat([1:(model.k2/data.Y)],model.N,1);
    ind23 = ind21+ind22+model.k1;
    ind24 = (ind23-1)*model.N + repmat([1:model.N]',1,size(ind23,2));
    ind2(ind24(:)) = 0;
    ind3 = find(ind2==2);
    alphaterm2(ind3) = 0; % zero out indices which are not active among latent topics
    indzero = [indzero; ind3];
end

if(model.phase==1 && model.option==6)
    ind2  = ones(model.N,model.k1);
    ind2  = [ind2 2*ones(model.N,model.K-model.k1)];
    ind21 = repmat(data.classlabels-1,1,model.k2)*model.k2;
    ind22 = repmat([1:model.k2],model.N,1);
    ind23 = ind21+ind22+model.k1;
    ind24 = (ind23-1)*model.N + repmat([1:model.N]',1,size(ind23,2));
    ind2(ind24(:)) = 0;
    ind3 = find(ind2==2);
    alphaterm2(ind3) = 0; % zero out indices which are not active among latent topics
    indzero = [indzero; ind3];
end

mod_alpha  = (M.*alphaterm1)*x + (1-x)*alphaterm2;   %alpha_{ni}  -- we need to ignore co-ordinates which are zero

gmod_alpha = gammaln(mod_alpha);
gmod_alpha(find(gmod_alpha==inf)) = 0;               %gammaln(alpha_{ni}) inf replaced by zero as we need to sum over the co-ordinates
der_alpha  = (M.*alphaterm1) - alphaterm2;

psi_alpha = psi(mod_alpha);
psi_alpha(find(psi_alpha==-inf)) = 0;

psi_gamma  = psi(model.gamma);
psi_gamma(find(psi_gamma==-inf)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term1 = sum(gammaln(sum(mod_alpha,2)));
term2 = - sum(sum(gmod_alpha));
term3 = (psi_gamma -repmat(psi(sum(model.gamma,2)),1,model.K)).*(mod_alpha-1);
term3(indzero) = 0;
f     = -(term1+term2+sum(sum(term3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    
    term11 = sum(psi(1,sum(mod_alpha,2)).*sum(der_alpha,2));
    term12 = -sum(sum(psi_alpha.*der_alpha));
    term13 = (psi_gamma -repmat(psi(sum(model.gamma,2)),1,model.K)).*(M.*alphaterm1-alphaterm2);
    term13(indzero) = 0;
    
    g = -(term11+term12+sum(sum(term13)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% function hessian %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % not correct; need to modify at some point
    if nargout > 2
        
        term21 = sum(psi(2,sum(mod_alpha,2)).*sum(der_alpha.^2,2));
        term22 = --sum(sum(psi(2,mod_alpha).*der_alpha));
        h = -diag(term21+term22);
        
    end
    
end

end
