function [model] = M_step(model, data, MAXMSTEPITER, MaxFun, maxvalue)

disp('M step starts');

sumphi = sum_phi(model, data); %% used in update of both epsilon and eta;

%% M step of DSLDA
if(model.option>=4)
    model   = update_alpha1_alpha2(model, data, MaxFun);
% %         disp('alpha1 done');
% %             value11M = cal_likelihood(model, data)
% %     
% %             if (compareval(value11M, maxvalue))
% %                 maxvalue = value11M;
% %             else
% %                 error('Incorrect after alpha1');
% %             end
    
    if(model.option~=7 && model.option~=8) %% exclude if DSLDA-OSST or DSLDA-NSLT2
        
        %% model   = update_epsilon(model, sumphi);
        % %         disp('epsilon done');
        % %                 value12M = cal_likelihood(model, data)
        % %
        % %                 if (compareval(value12M, maxvalue))
        % %                     maxvalue = value12M;
        % %                 else
        % %                     error('Incorrect after epsilon');
        % %                 end
        
    end
    
    
else %% for option = 1, 2 and 3
    model   = update_alpha(model, data, MaxFun);
    % %         value1M = cal_likelihood(model, data)
    % %
    % %         if (compareval(value1M, maxvalue))
    % %             maxvalue = value1M;
    % %         else
    % %             error('Incorrect after alpha');
    % %         end
end

model.log_beta = update_beta_cpp (model);
% % value2M     = cal_likelihood(model, data)
% % 
% % if (compareval(value2M, maxvalue))
% %     maxvalue = value2M;
% % else
% %     error('Incorrect after beta');
% % end

if(model.option>=3)
    model  = update_eta(model, data, sumphi);
    %%   value3 = cal_likelihood(model, data)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% don't include the following unless you calculate the slack variables from svm package and include that in the lower bound calculation
    %    maxvalue = value3;
    %     if (compareval(value3, maxvalue))
    %         maxvalue = value3;
    %     else
    %         error('Incorrect after eta');
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

disp('M step done');

end
