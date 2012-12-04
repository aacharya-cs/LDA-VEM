function [model] = update_alpha1 (model, data, MaxFun)

% if(model.option>=4)
%     options = optimset('LargeScale', 'on', 'Algorithm', 'interior-point', 'GradObj', 'on', 'MaxFunEvals', MaxFun, 'Display', 'off');
% else

% % options = optimset('Algorithm','trust-region-reflective','GradObj','on', 'MaxFunEvals', MaxFun, 'MaxIter', MaxFun, 'Display', 'off');
% % %end
% % 
% % x0      = [model.alpha1 zeros(1, model.K-model.k1)];
% % D       = length(x0);
% % lb      = model.MINVALUE*ones(1,D);
% % tmpvect = ones(1,D);
% % tmpvect(model.k1+1:end) = 0;
% % Aeq     = diag(ones(1,D));
% % beq     = zeros(1,D);
% % 
% % temp    = fmincon(@L_alpha1, x0, [], [], [], [], lb, [], [], options, model, data);
% % 
% % 
% % model.alpha1 = temp(1:model.k1);
% % tempver = model.alpha1 - x0(1:model.k1);
% % 
% % % % temp = dir_newton_alpha1(model.gamma, MaxFun, model.alpha);
% % % % model.alpha1 = temp(1:model.k1);

if(model.phase==1)
    M = [model.epsilon*data.annotations (1-model.epsilon)*ones(model.N, model.K-model.k1)];
else
    M = [model.epsilon*ones(model.N, model.k1) (1-model.epsilon)*ones(model.N,model.K-model.k1)];
end

if(model.phase==1 && model.option==5)
    ind2  = ones(model.N,model.k1);
    ind2  = [ind2 2*ones(model.N,model.K-model.k1)];
    ind21 = repmat(data.classlabels-1,1,(model.k2/data.Y))*(model.k2/data.Y);
    ind22 = repmat([1:(model.k2/data.Y)],model.N,1);
    ind23 = ind21+ind22+model.k1;
    ind24 = (ind23-1)*model.N + repmat([1:model.N]',1,size(ind23,2));
    ind2(ind24(:)) = 0;
    ind3 = find(ind2==2);
    M(ind3) = 0; % zero out indices which are not active among latent topics
end

if(model.phase==1 && model.option==6)
    ind2  = ones(model.N,model.k1);
    ind2  = [ind2 2*ones(model.N,model.K-model.k1)];
    ind21 = repmat(data.classlabels-1,1,model.k2)*model.k2;
    ind22 = repmat([1:model.k2],model.N,1);
    ind23 = ind21+ind22+model.k1;
    ind24 = (ind23-1)*model.N + repmat([1:model.N]',1,size(ind23,2));
    ind2(ind24(:)) = 0;
    ind3 = find(ind2==2);
    M(ind3) = 0; % zero out indices which are not active among latent topics
end

x0   = [model.alpha1 zeros(1, model.K-model.k1)];
temp = generalized_dir_newton_alpha (model.gamma, MaxFun, x0, M);
model.alpha1 = temp(1:model.k1);

end
