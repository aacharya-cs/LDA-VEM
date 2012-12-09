function [model] = update_lambda_m(model, option)

    temp = repmat(model.eta1, model.K1, 1);
    model.lambda(1:model.K1,:) = temp + model.ss_lambda(1:model.K1,:);
    if(option==1) %% NPDSLDA
      temp = repmat(model.eta2, model.K2, 1);
      model.lambda(model.K1+1:end,:) = temp + model.ss_lambda(model.K1+1:end,:);
    end

    %% normalize lambda
    model.lambda = model.lambda./repmat(sum(model.lambda,2),1,model.V);

end







