function [model] = optimal_ordering(model)

temp1 = model.ss_u;
temp2 = [1:model.K1];
temp3 = [temp1; temp2]';

%% sort by first column
[~,ind] = sort(temp3(:,1), 'descend');

%% ind = load('/lusr/u/ayan/Documents/DSLDA_SDM/DSLDA/onlineHDP/idx.txt');
%% ind = ind + 1;

temp3   = temp3(ind,:);
model.ss_u = temp3(:,1)';
model.ss_lambda(1:model.K1,:) = model.ss_lambda(ind,:);

end

