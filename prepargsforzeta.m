function [arg1 arg2 arg3] = prepargsforzeta(model)

term1 = psi(model.a);
term2 = cumsum(psi(model.b),2)  - psi(model.b);
term3 = cumsum(psi(model.a + model.b),2);
term4 = cumsum(psi(model.b),2) - cumsum(psi(model.a + model.b),2);
arg1  = (term1 + term2 - term3);
arg1(:,model.T) = term4(:,end); %% arg1 = E_{q}[log(\beta))

term5 = psi(model.lambda); %% psi(lambda_kv)
term6 = repmat(psi(sum(model.lambda,2)),1,model.V); %% psi(\sum_{v=1}^{V}lambda_kv)
arg2  = (term5 - term6);   %% psi(lambda_kv) - psi(\sum_{v=1}^{V}lambda_kv)

term7 = psi(model.mun);    %% psi(mu_nk)
ind   = find(term7==-inf);
term7(ind) = 0;
term8 = repmat(psi(sum(model.mun,2)),1,model.K2); %% psi(\sum_{k=1}^{K2}mu_nk)
arg3  = (term7 - term8);   %% psi(mu_nk) - psi(\sum_{k=1}^{K2}mu_nk)

end