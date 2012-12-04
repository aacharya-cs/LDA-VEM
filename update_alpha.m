function model = update_alpha (model, data, MaxFun)

% % options = optimset('Algorithm','trust-region-reflective','GradObj','on', 'MaxFunEvals', MaxFun, 'Display', 'off');
% % 
% % x0      = model.alpha;
% % D       = length(x0);
% % lb      = model.MINVALUE*ones(1,D);
% % temp    = fmincon(@L_alpha, x0, [], [], [], [], lb, [], [], options, model, data);
% % model.alpha = temp;

%% for Minka's method
%%model.alpha = dir_newton_alpha(model.gamma, MaxFun, model.alpha);



M    = [ones(model.N, model.K)];
x0   = model.alpha;

%temp2 = dir_newton_alpha (model.gamma, MaxFun, x0);
temp = generalized_dir_newton_alpha (model.gamma, MaxFun, x0, M);
model.alpha = temp;

end
