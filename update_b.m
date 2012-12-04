function [model] =  update_b(model)

temp1   = model.sumzeta;
temp1   = temp1(:,2:model.T);
temp2   = partial_sum(temp1);
model.b = model.alphazero(2) + temp2;

end