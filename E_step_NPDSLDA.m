function [model] = E_step_NPDSLDA(model, data, MAXESTEPITER, maxvalue, countVEM, option)

%% E step of proposed model

count1     = 0;
disp('E step starts');
%% maxvalue = likelihood_NPDSLDA(model, data, option)
%% this calculation is essential while debugging as the svm slack variable part is not included in lower bound calculation

while(count1< MAXESTEPITER)
    maxvalue = -1e50;
    disp('count from E-step');
    count2 = 0;
    while(count2<10)
        
        %%%%%%%%%% small-phi update  %%%%%%%%%%%%%%%%
        [arg1, arg2]   = prepargsforsmallphi(model);
        model.smallphi = update_smallphi(model, data, arg1, arg2, count2);
%         disp('small_phi update done');
%         value = likelihood_NPDSLDA(model, data, option)
%         if (compareval(value, maxvalue))
%             maxvalue = value;
%         else
%             error('Incorrect after small_phi');
%         end
        
        %%%%%%%%%% zeta update  %%%%%%%%%%%%%%%%
        term1 = squeeze(model.smallphi(1,:,:));
        [arg1, arg2, arg3] = prepargsforzeta(model);
        [model.zeta, model.ss_lambda, model.ss_features] = update_zeta(model, data, arg1, arg2, arg3, count2, option);
        model.sumzeta = sum_zeta(model, data, option);
%         disp('zeta update done');
%         value = likelihood_NPDSLDA(model, data, option)
%         if (compareval(value, maxvalue))
%             maxvalue = value;
%         else
%             error('Incorrect after zeta');
%         end
        
        %%%%%%%%%% a update  %%%%%%%%%%%%%%%%
        
        model = update_a(model, data);
%         disp('a update done');
%         value = likelihood_NPDSLDA(model, data, option)
%         if (compareval(value, maxvalue))
%             maxvalue = value;
%         else
%             error('Incorrect after a');
%         end
        
        %%%%%%%%%% b update  %%%%%%%%%%%%%%%%
        model = update_b(model);
%         disp('b update done');
%         value = likelihood_NPDSLDA(model, data, option)
%         if (compareval(value, maxvalue))
%             maxvalue = value;
%         else
%             error('Incorrect after b');
%         end
        count2 = count2+1;
    end

    %value = likelihood_NPDSLDA(model, data, option)
    %%%%%%%%%% update and re-order sufficient statistics for better performance %%%%%%%%%%%%%%%%
    
    model.ss_u = (squeeze(sum(sum(model.smallphi,1),2)))'; %% \sum_{j}\sum_{t}\phi_{jtk1}
    model = optimal_ordering(model);
    
    %%%%%%%%%% u update  %%%%%%%%%%%%%%%%
    model = update_u(model);
%     disp('u update done');
%     value = likelihood_NPDSLDA(model, data, option)
%     if (compareval(value, maxvalue))
%         maxvalue = value;
%     else
%         error('Incorrect after u');
%     end
    %%%%%%%%%% v update  %%%%%%%%%%%%%%%%
    
    model = update_v(model);
%     disp('v update done');
%     value = likelihood_NPDSLDA(model, data, option)
%     if (compareval(value, maxvalue))
%         maxvalue = value;
%     else
%         error('Incorrect after v');
%     end
    %%%%%%%%%% lambda update  %%%%%%%%%%%%%%%%
    model  = update_lambda_m(model, option);
%     disp('lambda update done');
%     value = likelihood_NPDSLDA(model, data, option)
%     if (compareval(value, maxvalue))
%         maxvalue = value;
%     else
%         error('Incorrect after lambda');
%     end
    
    %%%%%%%%%% mun update  %%%%%%%%%%%%%%%%
    if(option==1) %% NPDSLDA
        model  = update_mun(model, data);
        %disp('mu_n update done');
        %value = likelihood_NPDSLDA(model, data, option)
        %if (compareval(value, maxvalue))
        %    maxvalue = value;
        %else
        %    error('Incorrect after mu_n');
        %end
    end
    
    count1 = count1+1;
end

end

