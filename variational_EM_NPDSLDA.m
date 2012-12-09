function [model, perfval] = variational_EM_NPDSLDA (data, MAXCOUNT, MAXESTEPITER, MAXMSTEPITER, MaxFun,  p, K1, T, option, svmcval)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main code for running variational EM on NPDSLDA model
% @ Ayan Acharya, Date: 11.07.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[model, data] = init_params_NPDSLDA(data, K1, T, option, svmcval);
countVEM = 1;
maxvalue = -1e50;

while (countVEM<=MAXCOUNT)
    
    disp('count from V-EM');
    countVEM
    
    % %     if(countVEM>1 && model.option==3)
    % %         maxvalue   = likelihood_NPDSLDA(model, data, option);
    % %     end
    
    %% E step
    model    = E_step_NPDSLDA(model,data, MAXESTEPITER, maxvalue, countVEM, option);
%     value    = likelihood_NPDSLDA(model, data, option);
%     
%     %%for checking if lower bound increases after each iteration; useful for debugging -- should be commented out to save time while running experiments
%     if (compareval(value, maxvalue))
%         maxvalue = value;
%     else
%         error('Incorrect after E step');
%     end
    
    %% M step
    model    = M_step_NPDSLDA(model,data, MaxFun, maxvalue, option);
%     value    = likelihood_NPDSLDA(model, data, option);
%     
%     if (compareval(value, maxvalue))
%         maxvalue = value;
%     else
%         error('Incorrect after M step');
%     end
    countVEM = countVEM+1;
    
end


% calculate the accuracy on the training set
[accval, confmat] = cal_accuracy_NPDSLDA(model, data, option);

perfval.accval  = accval;
perfval.confmat = confmat;
perfval.multiclassacc = 100*sum(diag(confmat))/sum(sum(confmat));

end
