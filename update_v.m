function [model] =  update_v(model)

temp1   = model.ss_u(2:end);
temp1   = partial_sum(temp1);
model.v = model.gammazero(2) + temp1;

end