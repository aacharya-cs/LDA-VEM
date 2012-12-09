function [arg1] = Elog_globalsticks(model)

term1 = psi(model.u);
term2 = cumsum(psi(model.v))  - psi(model.v);
term3 = cumsum(psi(model.u + model.v));
term4 = cumsum(psi(model.v)) - cumsum(psi(model.u + model.v));
arg1  = (term1 + term2 - term3);
arg1(1,length(arg1)+1) = term4(end); %% arg1 = E_{q}[log(\beta))

end