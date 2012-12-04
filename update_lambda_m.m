function [model] = update_lambda_m(model)

temp = repmat(model.eta,(model.K1+model.K2),1);
model.lambda = temp + model.ss_lambda;
%% normalize lambda
model.lambda = model.lambda./repmat(sum(model.lambda,2),1,model.V);

end







