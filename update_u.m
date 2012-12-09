function [model] =  update_u(model)

temp1   = model.ss_u(1:end-1);
model.u = model.gammazero(1) + temp1;

end