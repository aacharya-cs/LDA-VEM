function [arg1, arg2] = prepargsforsmallphi(model)

term1 = psi(model.u);
term2 = cumsum(psi(model.v))  - psi(model.v);
term3 = cumsum(psi(model.u + model.v));
term4 = cumsum(psi(model.v)) - cumsum(psi(model.u + model.v));
arg1  = (term1 + term2 - term3);
arg1(1,model.K1) = term4(end); %% arg1 = E_{q}[log(\beta))

term5 = model.lambda(1:model.K1,:);
term6 = psi(term5);        %% psi(lambda_{kv})
term7 = repmat(psi(sum(term5,2)),1,model.V); %% psi(\sum_{v=1}^{V}lambda_{kv})
arg2  = (term6-term7);       %% E_{log(\eta)} = psi(lambda_{kv}) - psi(\sum_{v=1}^{V}lambda_{kv}) 

end