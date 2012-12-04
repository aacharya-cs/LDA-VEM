function model=update_beta (model, data)


for i=1:model.K
    temp = zeros(1,model.V);
    for j=1:model.V
        for n=1:model.N
            for m=1:data.nwordspdoc(n)
                if(data.w{n}(m)==j)
                    temp(j) = temp(j) + model.phi{n}(m,i);
                end
            end
        end       
    end
    model.beta(i,:) =  temp/sum(temp);
end

end
