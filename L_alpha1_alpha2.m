function [f, g, h] = L_alpha1_alpha2 (x, gammas, M)

[N, K] = size(M);
mod_alpha  = (M.*repmat(x,N,1));
gmod_alpha = gammaln(mod_alpha);
gmod_alpha(find(gmod_alpha==inf)) = 0;
psi_mod_alpha  = psi(mod_alpha);
psi_mod_alpha(find(psi_mod_alpha==-inf)) = 0;
psi_gamma  = psi(gammas);
psi_gamma(find(psi_gamma==-inf)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term1 = sum(gammaln(sum(mod_alpha,2))) - sum(sum(gmod_alpha));
term2 = (psi_gamma -repmat(psi(sum(gammas,2)),1,K)).*(mod_alpha-1);
ind = find(M==0);
term2(ind) = 0;
f     = -(term1+sum(sum(term2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    
    gterm1 = sum(M.*repmat(psi(sum(mod_alpha,2)),1,K),1);
    gterm2 = -sum(M.*psi_mod_alpha,1);
    gterm3 = sum((M.*[psi_gamma - repmat(psi(sum(gammas,2)),1,K)]),1);
    g = -(gterm1 + gterm2 + gterm3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% function hessian %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to do; not yet completed
    if nargout > 2
        
    end
    
end

end
