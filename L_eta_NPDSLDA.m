function [f, g] = L_eta_NPDSLDA(x0, model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term1 = (model.K1+model.K2)*(gammaln(sum(x0)) - sum(gammaln(x0)));
term2 = (x0-1).*(sum(psi(model.lambda),1) - repmat(sum(psi(sum(model.lambda,2))),1, model.V));
term3 = sum(term2);
f     = -(term1+term3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    
    gterm1 = (model.K1+model.K2)*(repmat(psi(sum(x0)),1,model.V) - psi(x0));
    gterm2 = (sum(psi(model.lambda),1) - sum(psi(sum(model.lambda,2))));
    g = -(gterm1 + gterm2);
    
end

end
