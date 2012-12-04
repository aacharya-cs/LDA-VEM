function [model] = sum_zeta_m(model, data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for n=1:model.N
    temp = model.zeta{n};
    temp = data.wcount{n}*temp;
    model.sumzeta(n,:) = temp;
end

end

