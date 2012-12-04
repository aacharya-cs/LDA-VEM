function [model] = update_smallphi_m(model, data, arg1, arg2)

Esticks       = arg1; %% Eq[log(\beta)]
Elogbeta      = arg2; %% psi(\lambda) - psi(\sum_{v=1}^{V} \lambda)

for i=1:model.N               %% loop over documents
    
    zetai   = model.zeta{i};  %% zeta_{i} double array
    windexi = data.windex{i}; %% pointer to windex{i}
    wcounti = data.wcount{i}; %% pointer to wcount{i}
    
    ndistWords = length(wcounti);
    
    for t=1:model.T            %% loop over topics (lower level truncation)
        logsum = 0;
        tmp = zeros(1,model.K1);
        for k1=1:model.K1      %% loop over topics (higher level truncation)
            value = 0;
            for w=1:ndistWords %% loop over (distinct) words
                value   = value + wcounti(w)*zetai(w,t)*Elogbeta(k1,windexi(w));
            end
            %%tmp(k1) = Esticks(k1);
            tmp(k1) = tmp(k1) + value;
            logsum  = log_sum(tmp(k1),logsum);
        end
        %% conversion from log space to real number
        for k1=1:model.K1
            if(logsum - tmp(k1)>50)
                model.smallphi(i,t,k1) = model.MINVALUE;
            end
            if(logsum - tmp(k1)<50)
                model.smallphi(i,t,k1) = exp(tmp(k1)-logsum)+model.MINVALUE;
            end
        end
        
    end
end

end

