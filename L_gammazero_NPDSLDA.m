function [f, g] = L_gammazero_NPDSLDA(x0, model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term1 = (model.K1-1)*(gammaln(sum(x0)) - sum(gammaln(x0)));
term2 = (x0(1)-1)*sum((psi(model.u) - psi(model.u+model.v)));
term3 = (x0(2)-1)*sum((psi(model.v) - psi(model.u+model.v)));
f     = -(term1+term2+term3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    
    gterm1 = (model.K1-1)*(repmat(psi(sum(x0)),1,2) - psi(x0));
    gterm2 = sum((psi(model.u) - psi(model.u+model.v)));
    gterm3 = sum((psi(model.v) - psi(model.u+model.v)));
    gterm4 = [gterm2 gterm3];
    g = -(gterm1 + gterm4);
    
end

end