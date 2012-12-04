function [model] = update_uzero_NPDSLDA(model, MaxFun)

options = optimset('LargeScale', 'on', 'Algorithm', 'interior-point', 'GradObj', 'on', 'MaxFunEvals', MaxFun, 'Display', 'off');
x0      = model.uzero;
D       = length(x0);
lb      = model.MINVALUE*ones(1,D);
temp    = fmincon(@L_uzero_NPDSLDA, x0, [], [], [], [], lb, [], [], options, model);
model.uzero = temp;

end