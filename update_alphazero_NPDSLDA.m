function [model] = update_alphazero_NPDSLDA(model, MaxFun)

options = optimset('LargeScale', 'on', 'Algorithm', 'interior-point', 'GradObj', 'on', 'MaxFunEvals', MaxFun, 'Display', 'off');
x0      = model.alphazero;
D       = length(x0);
lb      = model.MINVALUE*ones(1,D);
temp    = fmincon(@L_alphazero_NPDSLDA, x0, [], [], [], [], lb, [], [], options, model);
model.alphazero = temp;

end