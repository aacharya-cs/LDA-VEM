function [value] = likelihood_NPDSLDA(model, data, option)

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
%% initializing different terms for lower bound calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
term1p  = 0;
term2p  = 0;
term3p  = 0;
term4p  = 0;
term5p  = 0;
term6p  = 0;
term7p  = 0;
term8p  = 0;
term9p  = 0;
term10p = 0;
term11p = 0;
term12p = 0;
term13p = 0;
term14p = 0;
term15p = 0;
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
if(option==1)  %% for NPDSLDA
    term61 = sum(sum((model.sumzeta(:,1:model.T-1).*(log(1-model.epsilon) + psi(model.a) - psi(model.a+model.b)))));
    term62 = model.sumzeta(:,2:model.T);
    term62 = sum(sum((partial_sum(term62).*(psi(model.b) - psi(model.a + model.b)))));
    term63 = sum(sum((model.sumzeta(:,model.T+1:end).*(log(model.epsilon) + psimodmun - repmat(psi(sum(model.mun,2)),1,model.K2)))));
    term6p = term61+term62+term63;
else %% for NPLDA
    term61 = sum(sum((model.sumzeta(:,1:model.T-1).*(psi(model.a) - psi(model.a+model.b)))));
    term62 = model.sumzeta(:,2:model.T);
    term62 = sum(sum((partial_sum(term62).*(psi(model.b) - psi(model.a + model.b)))));
    term6p = term61+term62;
end


%% term8
if(option==1)  %% for NPDSLDA
    term81 = model.K2*(gammaln(sum(model.eta)) - sum(gammaln(model.eta)));
    term82 = model.lambda(model.K1+1:end,:);
    term83 = (model.eta-1).*(sum(psi(term82),1) - repmat(sum(psi(sum(term82,2))),1, model.V));
    term8p = (term81+sum(term83));
end

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
if(option==1)  %% for NPDSLDA
    term121 = model.lambda(model.K1+1:end,:);
    term122 = sum(gammaln(sum(term121,2))) - sum(sum(gammaln(term121)));
    term123 = sum(sum((term121-1).*(psi(term121) - repmat(psi(sum(term121,2)),1,model.V))));
    term12p = -(term122+term123);
end

%% term13
term13p = model.smallphi.*log(model.smallphi);
term13p = -sum(sum(sum(term13p)));

%% term14 and term7
term14p = 0;
term7p  = 0;
Eloglambda = [psi(model.lambda) - repmat(psi(sum(model.lambda,2)),1,model.V)];
for n=1:model.N
    %% term7
    term71 = model.zeta{n}(:,1:model.T)*squeeze(model.smallphi(n,:,:)); %% \sum_{t}\zeta_{nmt}*\varphi_{ntk1} ..dimension M_{n} x K_{1}
    if(option==1)
        term72 = model.zeta{n}(:,model.T+1:end);                            %% ..dimension M_{n} x K_{2}
    end
    term73 = repmat(model.r(data.classlabels(n),:),data.Y,1);
    Mn = length(data.windex{n});
    
    for w=1:Mn
        term711 = term71(w,:)*Eloglambda(1:model.K1,data.windex{n}(w)); %% \sum_{k_{1}}term71_{mk_{1}}*Eloglambda_{k_{1},w_{nm}} .. scalar
        term7p  = term7p + data.wcount{n}(w)*term711;
        if(option==1)
            term712 = term72(w,:)*Eloglambda(model.K1+1:end,data.windex{n}(w)); %% \sum_{k_{2}}term72_{mk_{2}}*Eloglambda_{k_{2},w_{nm}} .. scalar
            term7p  = term7p + data.wcount{n}(w)*term712;
        end
    end
    
    if(model.phase==1)
        term74  = (term73-model.r);                  %% r(y_{n}) - r(y) .. dimension Y x (K_{1}+K_{2})
        term74  = model.dmu(n,:)*term74;             %% dimension 1 x (K_{1}+K_{2})
        term7p  = term7p + (1/data.nwordspdoc(n))*sum(term74.*model.ss_features(n,:));
    end
    
    %% term14
    term141 = data.wcount{n};
    if(option==1) %% for NPDSLDA
        term142 = repmat(term141',1,(model.T+model.K2)).*model.zeta{n}.*log(model.zeta{n});
    else
        term142 = repmat(term141',1,model.T).*model.zeta{n}.*log(model.zeta{n});
    end
    ind     = find(isnan(term142)==1);
    term142(ind) = 0;
    term14p = term14p + sum(sum(term142));  
end
term14p = -term14p;

%% term15
if(option==1)
    term151 = sum(gammaln(sum(model.mun,2))) - sum(sum(gammalnmodmun));
    term152 = (model.mun-1).*(psimodmun - repmat(psi(sum(model.mun,2)),1,model.K2));
    term152(indzero) = 0;
    term152 = sum(sum(term152));
    term15p = -(term151+term152);
end

%% final value

value = term1p+term2p+term3p+term4p+term5p+term6p+term7p+term8p+term9p;
value = value+term10p+term11p+term12p+term13p+term14p+term15p;

if(isnan(value) || value==inf || value==-inf)
    keyboard;
end

end
