function [testmodel, perfval, probassign] = InferenceOnTest_NPDSLDA(trmodel, data, MAXESTEPITER)

%% inference on test data
probassign = [];
testmodel  = inittestmodel_NPDSLDA(trmodel, data);
count1     = 0;
maxvalue   = -1e25;

disp('inference on test data starts');

while(count1<20)
    
    disp('count from E-step in only inference');
        
        %%%%%%%%%% small-phi update  %%%%%%%%%%%%%%%%
        [arg1, arg2]   = prepargsforsmallphi(testmodel);
        testmodel.smallphi = update_smallphi(testmodel, data, arg1, arg2, count1);
        %disp('small_phi update done');
        %         value = likelihood_NPDSLDA(testmodel, data)
        %         if (compareval(value, maxvalue))
        %             maxvalue = value;
        %         else
        %             error('Incorrect after small_phi');
        %         end
        
        %%%%%%%%%% zeta update  %%%%%%%%%%%%%%%%
        
        [arg1, arg2, arg3] = prepargsforzeta(testmodel);
        [testmodel.zeta, testmodel.ss_lambda] = update_zeta(testmodel, data, arg1, arg2, arg3, count1);
        testmodel.sumzeta = sum_zeta(testmodel, data);
        %disp('zeta update done');
        %         value = likelihood_NPDSLDA(testmodel, data)
        %         if (compareval(value, maxvalue))
        %             maxvalue = value;
        %         else
        %             error('Incorrect after zeta');
        %         end
        
        %%%%%%%%%% a update  %%%%%%%%%%%%%%%%
        
        testmodel = update_a(testmodel, data);
        %disp('a update done');
        %         if(count1==0 && count2==0)
        %         else
        %             value = likelihood_NPDSLDA(testmodel, data)
        %             if (compareval(value, maxvalue))
        %                 maxvalue = value;
        %             else
        %                 error('Incorrect after a');
        %             end
        %         end
        %%%%%%%%%% b update  %%%%%%%%%%%%%%%%
        testmodel = update_b(testmodel);
        %disp('b update done');
        %         if(count1==0 && count2==0)
        %         else
        %             value = likelihood_NPDSLDA(testmodel, data)
        %             if (compareval(value, maxvalue))
        %                 maxvalue = value;
        %             else
        %                 error('Incorrect after b');
        %             end
        %         end
    
    testmodel  = update_mun(testmodel, data);
    %disp('mu_n update done');
    
    %     value = likelihood_NPDSLDA(testmodel, data)
    %     if (compareval(value, maxvalue))
    %         maxvalue = value;
    %     else
    %         error('Incorrect after mu_n');
    %     end
    count1 = count1+1;
    
end

[accval, confmat, probassign, binacc] = cal_accuracy_NPDSLDA(testmodel, data);

perfval.acc = accval;
[I,J,K] =  unique(data.classlabels,'first');
Jprime = [J(2:end); length(data.classlabels)];
clscount = (Jprime-J);
perfval.wacc = 100*(perfval.acc')*clscount/sum(clscount);
perfval.confmat = confmat;
perfval.multiclassacc = 100*sum(diag(confmat))/sum(sum(confmat));
perfval.binacc = binacc;

end

