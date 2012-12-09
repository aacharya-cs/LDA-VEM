function [alpha, lambda] = hdp_to_lda(trmodel)

%% alpha
sticks = trmodel.u./(trmodel.u + trmodel.v);
alpha  = zeros(trmodel.K1,1);
left   = 1.0;
for i=1:(trmodel.K1-1)
    alpha(i) = sticks(i) * left;
    left     = left - alpha(i);
end
alpha(trmodel.K1-1) = left;
alpha = alpha';
alpha = alpha*10;
%%alpha = alpha*model.alphazero(1)/model.alphazero(1);

%% lambda
lambda = trmodel.lambda;

end