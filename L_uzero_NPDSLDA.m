function [f, g] = L_uzero_NPDSLDA(x0, model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psimodmun  = psi(model.mun);
psimodmun(find(psimodmun==-inf)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term1 = (model.N)*[gammaln(sum(x0)) - sum(gammaln(x0))];
term2 = (x0-1).*(sum(psimodmun,1) - repmat(sum(psi(sum(model.mun,2))),1,model.K2));
term3 = sum(term2);
f     = -(term1+term3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    
    gterm1 = (model.N)*(repmat(psi(sum(x0)),1,model.K2) - psi(x0));
    gterm2 = (sum(psimodmun,1) - sum(psi(sum(model.mun,2))));
    g = -(gterm1 + gterm2);
    
end

end
