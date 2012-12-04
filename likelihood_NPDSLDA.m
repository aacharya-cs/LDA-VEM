function [value] = likelihood_NPDSLDA(model, data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psimodmun  = psi(model.mun);
psimodmun(find(psimodmun==-inf)) = 0;

gammalnmodmun = gammaln(model.mun);
gammalnmodmun(find(gammalnmodmun==inf)) = 0;

mod_uzero  = (data.annotations.*repmat(model.uzero, model.N, 1));
gammaln_mod_uzero = gammaln(mod_uzero);
gammaln_mod_uzero(find(gammaln_mod_uzero==inf)) = 0;

indzero = find(data.annotations==0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% term1
term11 = model.K1*(gammaln(sum(model.eta)) - sum(gammaln(model.eta)));
term12 = model.lambda(1:model.K1,:);
term13 = (model.eta-1).*(sum(psi(term12),1) - repmat(sum(psi(sum(term12,2))),1, model.V));
term14 = sum(term13);
term1p = (term11+term14);

%% term2
term21 = sum(gammaln(sum(mod_uzero,2))) - sum(sum(gammaln_mod_uzero));
term22 = (mod_uzero-1).*(psimodmun - repmat(psi(sum(model.mun,2)),1,model.K2));
term22(indzero) = 0;
term2p = (term21+sum(sum(term22)));

%% term3
term31 = (model.K1-1)*(gammaln(sum(model.gammazero)) - sum(gammaln(model.gammazero)));
term32 = (model.gammazero(1)-1)*(psi(model.u) - psi(model.u + model.v));
term33 = (model.gammazero(2)-1)*(psi(model.v) - psi(model.u + model.v));
term3p = term31+sum(term32+term33);

%% term4
term41 = (model.N*(model.T-1))*(gammaln(sum(model.alphazero)) - sum(gammaln(model.alphazero)));
term42 = (model.alphazero(1)-1)*(psi(model.a) - psi(model.a+model.b));
term43 = (model.alphazero(2)-1)*(psi(model.b) - psi(model.a+model.b));
term4p = term41+sum(sum(term42+term43));

%% term5
term51 = (squeeze(sum(sum(model.smallphi,1),2)))';
term52 = partial_sum(term51(2:end)).*(psi(model.v) - psi(model.u + model.v));
term53 = (squeeze(sum(sum(model.smallphi,1),2)))';
term54 = term53(1:end-1).*(psi(model.u) - psi(model.u + model.v));
term5p = sum(term52+term54);

%% term6
term61 = sum(sum((model.sumzeta(:,1:model.T-1).*(log(1-model.epsilon) + psi(model.a) - psi(model.a+model.b)))));
term62 = model.sumzeta(:,2:model.T);
term62 = sum(sum((partial_sum(term62).*(psi(model.b) - psi(model.a + model.b)))));
term63 = sum(sum((model.sumzeta(:,model.T+1:end).*(log(model.epsilon) + psimodmun - repmat(psi(sum(model.mun,2)),1,model.K2)))));
term6p = term61+term62+term63;

%% term8
term81 = model.K2*(gammaln(sum(model.eta)) - sum(gammaln(model.eta)));
term82 = model.lambda(model.K1+1:end,:);
if(min(min(term82))<0)
    ind = find(term82(model.K1+1:end,:)<0);
    keyboard;
end
term83 = (model.eta-1).*(sum(psi(term82),1) - repmat(sum(psi(sum(term82,2))),1, model.V));
term84 = sum(term83);
term8p = (term81+term84);

%% term9
term91 = gammaln(model.u+model.v) - gammaln(model.u) - gammaln(model.v);
term92 = (model.u-1).*(psi(model.u)-psi(model.u+model.v));
term93 = (model.v-1).*(psi(model.v)-psi(model.u+model.v));
term9p = -sum(term91+term92+term93);

%% term10
term101 = gammaln(model.a+model.b) - gammaln(model.a) - gammaln(model.b);
term102 = (model.a-1).*(psi(model.a)-psi(model.a+model.b));
term103 = (model.b-1).*(psi(model.b)-psi(model.a+model.b));
term10p = -sum(sum(term101+term102+term103));

%% term11
term111 = model.lambda(1:model.K1,:);
term112 = sum(gammaln(sum(term111,2))) - sum(sum(gammaln(term111)));
term113 = sum(sum((term111-1).*(psi(term111) - repmat(psi(sum(term111,2)),1,model.V))));
term11p = -(term112+term113);

%% term12
term121 = model.lambda(model.K1+1:end,:);
term122 = sum(gammaln(sum(term121,2))) - sum(sum(gammaln(term121)));
term123 = sum(sum((term121-1).*(psi(term121) - repmat(psi(sum(term121,2)),1,model.V))));
term12p = -(term122+term123);

%% term13
term13p = model.smallphi.*log(model.smallphi);
ind     = find(isnan(term13p)==1);
term13p(ind) = 0;
term13p = -sum(sum(sum(term13p)));

%% term14 and term7
term14p = 0;
term7p  = 0;
for n=1:model.N    
    %% term7
    term71 = model.zeta{n}(:,1:model.T)*squeeze(model.smallphi(n,:,:));
    term72 = model.zeta{n}(:,model.T+1:end);
    term73 = repmat(model.r(data.classlabels(n),:),data.Y,1);
    for w=1:length(data.windex(n))
        term7p = term7p + data.wcount{n}(w)*(term71(w,:)*(psi(model.lambda(1:model.K1,data.windex{n}(w))) - psi(sum(model.lambda(1:model.K1,:),2))));
        term7p = term7p + data.wcount{n}(w)*(term72(w,:)*(psi(model.lambda(model.K1+1:end,data.windex{n}(w))) - psi(sum(model.lambda(model.K1+1:end,:),2))));
    end
    if(model.phase==1)
        term74  = (term73-model.r);
        term74  = model.dmu(n,:)*term74;
        term7p  = term7p + (1/data.nwordspdoc(n))*sum(term74.*model.sumzeta(n,:));
    end
    
    %% term14
    term141 = model.zeta{n}.*log(model.zeta{n}); 
    ind     = find(isnan(term141)==1); 
    term141(ind) = 0;
    term14p = term14p + sum(sum(term141)); 
    
end
term14p = -term14p;

%% term15
term151 = sum(gammaln(sum(model.mun,2))) - sum(sum(gammalnmodmun));
term152 = (model.mun-1).*(psimodmun - repmat(psi(sum(model.mun,2)),1,model.K2));
term152(indzero) = 0;
term152 = sum(sum(term152));
term15p = -(term151+term152);

%% final value

value = term1p+term2p+term3p+term4p+term5p+term6p+term7p+term8p+term9p;
value = value+term10p+term11p+term12p+term13p+term14p+term15p;

if(isnan(value) || value==inf || value==-inf)
    keyboard;
end

end
