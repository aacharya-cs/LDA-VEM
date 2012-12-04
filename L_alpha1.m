function [f, g, h] = L_alpha1 (x, model, data)

indzero = [];
x = x(1:model.k1);

if(model.phase==1)
    M = [data.annotations ones(model.N, model.K-model.k1)];
else
    M = [ones(model.N, model.K)];
end

if(model.phase==1 && model.option>=4)
    indzero = find(M==0);
end

alphaterm1 = [repmat(x,model.N,1) zeros(model.N, model.K-model.k1)];
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

mod_alpha  = (M.*alphaterm1)*model.epsilon + (1-model.epsilon)*alphaterm2;   %alpha_{ni}  -- we need to ignore co-ordinates which are zero

gmod_alpha = gammaln(mod_alpha);
gmod_alpha(find(gmod_alpha==inf)) = 0;               %gammaln(alpha_{ni}) inf replaced by zero as we need to sum over the co-ordinates
der_alpha  = (M.*alphaterm1) - alphaterm2;

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

%% for debugging
valuetemp = cal_likelihood(model, data); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    
    term11 = sum(psi(sum(mod_alpha,2)));
    term12 = -sum(psi_mod_alpha,1);
    term13 = (psi_gamma -repmat(psi(sum(model.gamma,2)),1,model.K))*(1-model.epsilon);
    term13(indzero) = 0;
    g = -(term11+term12+sum(term13,1));
    g(model.k1+1:end) = 0;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% function hessian %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to do; not yet completed
    if nargout > 2
        
        term21 = sum(psi(1,sum(mod_alpha,2)));
        term22 = -sum(psi(1,mod_alpha),1);
        h = -diag(term21+term22)*(1-model.epsilon)^2;
        
    end
    
end

end
