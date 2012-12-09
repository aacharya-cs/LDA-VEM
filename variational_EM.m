function [model, perfval] = variational_EM (data, MAXCOUNT, MAXESTEPITER, MAXMSTEPITER, MaxFun,  p, K, option, phase, model, epsilon, svmcval)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main code for running variational EM on the proposed model
% @ Ayan Acharya, Date: 05.28.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model    = init_params(data, K, data.V, p, option, 1, epsilon, svmcval);
countVEM = 1;
maxvalue = -100000000000000000000;

while (countVEM<=MAXCOUNT)
    
    disp('count from V-EM');
    countVEM
    
    % %     if(countVEM>1 && model.option==3)
    % %         maxvalue   = cal_liKelihood(model, data);
    % %     end
    
    %% E step
    model    = E_step(model,data, MAXESTEPITER, maxvalue, option, countVEM);
    % %     value1   = cal_liKelihood(model, data)
    % %
    % %     %%for checking if lower bound increases after each iteration; useful for debugging -- should be commented out to save time while running experiments
    % %     if (compareval(value1, maxvalue))
    % %         maxvalue = value1;
    % %     else
    % %         error('Incorrect after E step');
    % %     end
    
    %% M step
    model    = M_step(model,data, MAXMSTEPITER, MaxFun, maxvalue);
    %% value2   = cal_liKelihood(model, data);
    
    %% if (compareval(value2, maxvalue))
    %%    maxvalue = value2;
    %% else
    %%    error('Incorrect after M step');
    %% end
    countVEM = countVEM+1;
    
end

% calculate the accuracy on the training set
[accval, confmat] = cal_accuracy(model, data);

perfval.accval  = accval;
perfval.confmat = confmat;
perfval.multiclassacc = sum(diag(confmat))/sum(sum(confmat));

end
