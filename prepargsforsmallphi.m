function [arg1, arg2] = prepargsforsmallphi(model)

arg1  = Elog_globalsticks(model);
term5 = model.lambda(1:model.K1,:);
term6 = psi(term5);          %% psi(lambda_{kv})
term7 = repmat(psi(sum(term5,2)),1,model.V); %% psi(\sum_{v=1}^{V}lambda_{kv})
arg2  = (term6-term7);       %% E_{log(\eta)} = psi(lambda_{kv}) - psi(\sum_{v=1}^{V}lambda_{kv}) 

end