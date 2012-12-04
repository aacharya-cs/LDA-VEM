function [model] = E_step(model, data, MAXESTEPITER, maxvalue, option, countVEM)

% E step of proposed model

count     = 0;
disp('E step starts');

%% this calculation is essential while debugging as the svm slack variable part is not included in lower bound calculation
if(countVEM>1)
 %%maxvalue = cal_likelihood(model, data);
end

if(countVEM>=1) %% have to change
while(count<MAXESTEPITER)
    
    %disp('count from E-step');
    count  = count+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% for debugging
    %%% model = update_phi(model, data);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (option==1)        %% unsupervised LDA
        [model.phi model.ss_topicword model.ss_topic] = update_phi_cpp(model, data, psi(model.gamma), option, [], []);
    elseif (option==2)    %% labeled LDA
        [model.phi model.ss_topicword model.ss_topic] = update_phi_cpp(model, data, psi(model.gamma), option, 1, data.annotations);
    elseif (option==3)    %% medLDA
        [model.phi model.ss_topicword model.ss_topic] = update_phi_cpp(model, data, psi(model.gamma), option, 1, []);
    elseif (option>=4)    %% DSLDA(4), DLSDA-NSLT-Ray(5), DSLDA-NSLT-Ayan(6), DSLDA-OSST(7)
        [model.phi model.ss_topicword model.ss_topic] = update_phi_cpp(model, data, psi(model.gamma), option, 1, [data.annotations ones(model.N,model.k2)]);
    else
    end
    
    %% for checking if lower bound increases after each update; useful for debugging -- should be commented out to save time while running experiments
% %     disp('phi update done');
% %     value1 = cal_likelihood(model, data)
% %     if (compareval(value1, maxvalue))
% %         maxvalue = value1;
% %     else
% %         error('Incorrect after phi');
% %     end
    
    model  = update_gamma(model, data);
% %     disp('gamma update done');    
% %     value2 = cal_likelihood(model, data)
% %     if (compareval(value2, maxvalue))
% %         maxvalue = value2;
% %     else
% %         error('Incorrect after gamma');
% %     end
end

else
 disp('skipping the first E-step .. EM will start from M step');
end

end

