function model=update_phi(model, data)

psi_gamma = psi(model.gamma);
psi_gamma(find(psi_gamma==-inf)) = 0;

for n=1:model.N
    val = (1/data.nwordspdoc(n))*[model.mu(n,:)*(repmat(model.eta(data.classlabels(n),:),model.Y,1) - model.eta)];
    for m=1:data.nwordspdoc(n)
        for i=1:model.K
            temp(i) = psi_gamma(n,i) + log(model.beta(i,data.w{n}(m))) - val(i);
        end
        logsum=0;
        for i=1:model.K
            logsum = sumlog(temp(i),logsum);
        end
        if(logsum==0)
            error('problem with losum=0');
        end
        for i=1:model.K
            if(logsum-temp(i)>50)
                model.phi{n}(m,i) = model.MINVALUE;
            else
                model.phi{n}(m,i) = exp(temp(i)-logsum) + model.MINVALUE;
            end
        end
        ind1 = isnan(model.phi{n}(m,:));
        if(~isempty(find(ind1==1)))
            error('problem with nan');
        end
    end
    ind2 = isnan(val);
    if(~isempty(find(ind2==1)))
        error('problem with val becoming nan');
    end
end

end

