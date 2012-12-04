function [model] = update_gammazero_NPDSLDA(model, MaxFun)

options = optimset('LargeScale', 'on', 'Algorithm', 'interior-point', 'GradObj', 'on', 'MaxFunEvals', MaxFun, 'Display', 'off');
x0      = model.gammazero;
D       = length(x0);
lb      = model.MINVALUE*ones(1,D);
temp    = fmincon(@L_gammazero_NPDSLDA, x0, [], [], [], [], lb, [], [], options, model);
model.gammazero = temp;

end