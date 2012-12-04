function [model] =  update_a(model, data)

temp1   = model.sumzeta;
temp1   = temp1(:,1:model.T-1);
model.a = model.alphazero(1) + temp1;

end