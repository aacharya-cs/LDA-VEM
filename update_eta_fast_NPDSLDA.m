function [eta] = update_eta_fast_NPDSLDA(lambda, MaxIter, ini_eta)

%% @Ayan Acharya
%% eta update using Minka's code

[M,K] = size(lambda);

if nargin < 3
    ini_eta = mean(lambda) / K; %% initial point
    if nargin < 2
        maxiter = 20;
    end
end

l = 0;
g = zeros(1,K);
pg = sum(psi(lambda),1) - sum(psi(sum(lambda,2)));
eta = ini_eta;
peta = zeros(1,K);

for t = 1:MaxIter
    l = l + 1;
    eta0 = sum(eta);
    g = M * (psi(eta0) - psi(eta)) + pg;
    h = - 1 ./ psi(1,eta);
    hgz = h * g' / (1 / psi(1,eta0) + sum(h));    
    eta = eta - h.*(g - hgz)/M;
    if any(eta < 0)
        eta = update_eta_fast_NPDSLDA(lambda, MaxIter, ini_eta/10); % try again!
        return;
    end  
    if (l > 1) && converged(eta,peta,1.0e-7)
        break;
    end
    peta = eta;
end


end
