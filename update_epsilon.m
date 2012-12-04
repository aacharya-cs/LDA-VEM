function [model] = update_epsilon (model, sumphi)

model.epsilon = sum(sum(sumphi(:,1:model.k1)))/sum(sum(sumphi));

end
