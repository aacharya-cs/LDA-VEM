function [model] = E_step_NPDSLDA(model, data, MAXESTEPITER, maxvalue, countVEM)

%% E step of proposed model

count1     = 0;
disp('E step starts');
maxvalue = -1e50; %%likelihood_NPDSLDA(model, data)
%% this calculation is essential while debugging as the svm slack variable part is not included in lower bound calculation

while(count1< MAXESTEPITER)
    maxvalue = -1e50;
    disp('count from E-step');
    count2 = 0;
    while(count2<4)
        
        %%%%%%%%%% small-phi update  %%%%%%%%%%%%%%%%
        [arg1, arg2]   = prepargsforsmallphi(model);
        model.smallphi = update_smallphi(model, data, arg1, arg2, count2);
        %         disp('small_phi update done');
        %         value = likelihood_NPDSLDA(model, data)
        %         if (compareval(value, maxvalue))
        %             maxvalue = value;
        %         else
        %             error('Incorrect after small_phi');
        %         end
        
        %%%%%%%%%% zeta update  %%%%%%%%%%%%%%%%
        
        [arg1, arg2, arg3] = prepargsforzeta(model);
        [model.zeta, model.ss_lambda] = update_zeta(model, data, arg1, arg2, arg3, count2);
        model.sumzeta = sum_zeta(model, data);
        %         disp('zeta update done');
        %         value = likelihood_NPDSLDA(model, data)
        %         if (compareval(value, maxvalue))
        %             maxvalue = value;
        %         else
        %             error('Incorrect after zeta');
        %         end
        
        %%%%%%%%%% a update  %%%%%%%%%%%%%%%%
        
        model = update_a(model, data);
        %disp('a update done');
        %         if(count1==0 && count2==0)
        %         else
        %             value = likelihood_NPDSLDA(model, data)
        %             if (compareval(value, maxvalue))
        %                 maxvalue = value;
        %             else
        %                 error('Incorrect after a');
        %             end
        %         end
        %%%%%%%%%% b update  %%%%%%%%%%%%%%%%
        model = update_b(model);
        %disp('b update done');
        %         if(count1==0 && count2==0)
        %         else
        %             value = likelihood_NPDSLDA(model, data)
        %             if (compareval(value, maxvalue))
        %                 maxvalue = value;
        %             else
        %                 error('Incorrect after b');
        %             end
        %         end
        count2 = count2+1;
    end
    %value = likelihood_NPDSLDA(model, data)
    %%%%%%%%%% update and re-order sufficient statistics for better performance %%%%%%%%%%%%%%%%
    
    %     temp = zeros(model.K1, model.V);
    %     for k=1:model.K1
    %         for n=1:model.N
    %             for w=1:length(data.windex{n})
    %                 tempvar = squeeze(model.smallphi(n,:,k))*(squeeze(model.zeta{n}(w,1:model.T)))';
    %                 temp(k,data.windex{n}(w)) = temp(k,data.windex{n}(w)) + data.wcount{n}(w)*tempvar;
    %             end
    %         end
    % %     end
    %
    %
    %     B = model.ss_lambda(1:model.K1,:);
    %     A = load('/home/ayan/Dropbox/DSLDA_SDM/onlineHDP/var_beta_ss.txt');
    %     C = max(max(A-B));
    %     D = max(max(A-temp));
    
    model.ss_u = (squeeze(sum(sum(model.smallphi,1),2)))'; %% \sum_{j}\sum_{t}\phi_{jtk1}
    model = optimal_ordering(model);
    
    %     BB = model.ss_lambda(1:model.K1,:);
    %     AA = load('/home/ayan/Dropbox/DSLDA_SDM/onlineHDP/var_beta_ss2.txt');
    %     CC = max(max(AA-BB));
    
    %%%%%%%%%% u update  %%%%%%%%%%%%%%%%
    model = update_u(model);
%         disp('u update done');
%         value = likelihood_NPDSLDA(model, data)
%         if (compareval(value, maxvalue))
%             maxvalue = value;
%         else
%             error('Incorrect after u');
%         end
    %%%%%%%%%% v update  %%%%%%%%%%%%%%%%
    
    model = update_v(model);
%         disp('v update done');
%         value = likelihood_NPDSLDA(model, data)
%         if (compareval(value, maxvalue))
%             maxvalue = value;
%         else
%             error('Incorrect after v');
%         end
    %%%%%%%%%% lambda update  %%%%%%%%%%%%%%%%
    model = update_lambda_m(model);
%         disp('lambda update done');
%         value = likelihood_NPDSLDA(model, data)
%         if (compareval(value, maxvalue))
%             maxvalue = value;
%         else
%             error('Incorrect after lambda');
%         end
    
    %%%%%%%%%% mun update  %%%%%%%%%%%%%%%%
    model  = update_mun(model, data);
%         disp('mu_n update done');
%         value = likelihood_NPDSLDA(model, data)
%         if (compareval(value, maxvalue))
%             maxvalue = value;
%         else
%             error('Incorrect after mu_n');
%         end
    
    count1 = count1+1;
end

end

