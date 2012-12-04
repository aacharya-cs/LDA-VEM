function [term3] = L_eta(model,data)

term3 = 0;
tempterm3 = diag(ones(1,k))*(-1)+1;
for l=1:model.r1
    tempw1 = data.dataw1(l).w;
    term3  = term3 - sum(sum(tempw1.*(tempterm3.model.eta*model.gamma')'));
end

term3 = -term3; % returns negative of the actual value

end
