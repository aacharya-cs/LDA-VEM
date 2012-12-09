function [alpha, lambda] = hdp_to_lda(model)

        %% alpha 
        sticks = (model.u)./(mode.u+model.v);
        alpha  = zeros(model.K1);
        left   = 1.0;
        for i:model.K1
            alpha[i] = sticks[i] * left;
            left = left - alpha[i];
        end

        alpha[model.K1-1] = left;      
        alpha = alpha*model.alphazero(1)/model.alphazero(2);
        
        %% lambda
        lambda_sum = sum(model.lambda, 2);
        lambda = model.lambda./repmat(lambda_sum,1,model.V);

end
