function [model] = update_eta_NPDSLDA(model, MaxFun)

options = optimset('LargeScale', 'on', 'Algorithm', 'interior-point', 'GradObj', 'on', 'MaxFunEvals', MaxFun, 'Display', 'off');
x0      = model.eta;
D       = length(x0);
lb      = model.MINVALUE*ones(1,D);
temp    = fmincon(@L_eta_NPDSLDA, x0, [], [], [], [], lb, [], [], options, model);
model.eta = temp;

end
