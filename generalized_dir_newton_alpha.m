function [alpha] = generalized_dir_newton_alpha (gammas, maxiter, ini_alpha, M)

[N, K] = size(gammas);
l = 0;
psi_gamma = psi(gammas);
psi_gamma(find(psi_gamma==-inf)) = 0;
gterm3 = sum((M.*[psi_gamma - repmat(psi(sum(gammas,2)),1,K)]),1);
alpha = ini_alpha;
palpha = zeros(1,K);
g = zeros(1,K);


for t = 1:maxiter
    
    l = l + 1;
 
    mod_alpha = M.*repmat(alpha,N,1);
    psi_mod_alpha = psi(mod_alpha);
    psi_mod_alpha(find(psi_mod_alpha==-inf)) = 0;
    tri_mod_alpha = psi(1,mod_alpha);
    tri_mod_alpha(find(tri_mod_alpha==inf))=0;

    gterm1 = sum(M.*repmat(psi(sum(mod_alpha,2)),1,K),1);    
    gterm2 = -sum(M.*psi_mod_alpha,1);
    
    g = gterm1 + gterm2 + gterm3;
    
    qinv = -1./sum((M.*M).*tri_mod_alpha,1);
    qinv(find(qinv==-inf)) = 0; %% not sure yet
    
    zinvterm1 = (M.*M);
    zinvterm2 = psi(1,repmat(sum(mod_alpha,2),1,K));    
    zinv = 1./sum(zinvterm1.*zinvterm2,1);
    zinv(find(zinv==inf)) = 0; %% not sure yet
    
    b = sum(g.*qinv)./(zinv+sum(qinv));
    
    hginv = (g-b).*qinv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% % for debugging    
% %     alpha0 = sum(alpha);
% %     g1 = N * (psi(alpha0) - psi(alpha)) + pg;
% %     h1 = - 1 ./ psi(1,alpha);
% %     hgz = h1 * g1' / (1 / psi(1,alpha0) + sum(h1));
% %     hginv1 = h1.*(g1 - hgz) / N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    alpha = alpha - hginv;
    
    if any(alpha < 0)
        alpha = generalized_dir_newton_alpha (gammas,maxiter,ini_alpha/10, M); % try again!
        return;
    end
    if (l > 1) && converged(alpha,palpha,1.0e-4)
        break;
    end
    palpha = alpha;
end

end

