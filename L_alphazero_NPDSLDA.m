function [f, g] = L_alphazero_NPDSLDA(x0, model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term1 = (model.N*(model.T-1))*(gammaln(sum(x0)) - sum(gammaln(x0)));
term2 = (x0(1)-1)*sum(sum((psi(model.a) - psi(model.a+model.b))));
term3 = (x0(2)-1)*sum(sum((psi(model.b) - psi(model.a+model.b))));
f     = -(term1+term2+term3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    
    gterm1 = (model.N*(model.T-1))*(repmat(psi(sum(x0)),1,2) - psi(x0));
    gterm2 = sum(sum((psi(model.a) - psi(model.a+model.b))));
    gterm3 = sum(sum((psi(model.b) - psi(model.a+model.b))));
    gterm4 = [gterm2 gterm3];
    g = -(gterm1 + gterm4);
    
end

end