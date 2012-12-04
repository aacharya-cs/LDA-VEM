function [testmodel, perfval, probassign] = InferenceOnTest(trmodel, data, MAXESTEPITER, option, p)

%% inference on test data
probassign = [];
testmodel  = inittestmodel(trmodel, data, option, p);
count      = 0;
maxvalue   = -100000000000000000000;

disp('inference on test data starts');

while(count<MAXESTEPITER)
    
    disp('count from E-step in only inference');
    count  = count+1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% for debugging 
    %%% model = update_phi(model, data);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if (option==1)        %% unsupervised LDA
        testmodel.phi = update_phi_cpp(testmodel, data, psi(testmodel.gamma), 1, [], []);
    elseif (option==2)    %% labeled LDA
        [testmodel.phi testmodel.ss_topicword testmodel.ss_topic] = update_phi_cpp(testmodel, data, psi(testmodel.gamma), 2, 0, data.annotations);
    elseif (option==3)    %% medLDA
        [testmodel.phi testmodel.ss_topicword testmodel.ss_topic] = update_phi_cpp(testmodel, data, psi(testmodel.gamma), 3, 0, []);
    elseif (option>=4)    %% DSLDA, DLSDA-NSLT-Ray, DSLDA-NSLT-Ayan, DSLDA-OSST
        [testmodel.phi testmodel.ss_topicword testmodel.ss_topic]= update_phi_cpp(testmodel, data, psi(testmodel.gamma), option, 0, [data.annotations ones(testmodel.N,testmodel.k2)]); %% to do
    else
    end
    
    %% for checking if lower bound increases after each update; useful for debugging -- should be commented out to save time while running experiments
% %     value1 = cal_likelihood(testmodel, data)
% %     if (compareval(value1, maxvalue))
% %        maxvalue = value1;
% %     else
% %        error('Incorrect after phi');
% %     end
    
    testmodel  = update_gamma(testmodel, data);
    
% %     value2 = cal_likelihood(testmodel, data) 
% %     if (compareval(value2, maxvalue))
% %        maxvalue = value2;
% %     else
% %        error('Incorrect after gamma');
% %     end
    
end

[accval, confmat, probassign, binacc] = cal_accuracy(testmodel, data);

perfval.acc = accval;
[I,J,K] =  unique(data.classlabels,'first');
Jprime = [J(2:end); length(data.classlabels)];
clscount = (Jprime-J);
perfval.wacc = 100*(perfval.acc')*clscount/sum(clscount);
perfval.confmat = confmat;
perfval.multiclassacc = 100*sum(diag(confmat))/sum(sum(confmat));
perfval.binacc = binacc;

end

