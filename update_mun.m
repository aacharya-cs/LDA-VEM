function [model] =  update_mun(model, data)

temp1 = model.sumzeta;
temp1 = temp1(:,model.T+1:end);
temp2 = repmat(model.uzero, model.N, 1) + temp1;

if(model.phase==1)
    temp3 = temp2.*data.annotations;
    model.mun = temp3;
end

end