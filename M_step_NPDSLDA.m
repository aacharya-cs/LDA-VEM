function [model] = M_step(model, data, MaxFun, maxvalue, option)

disp('M step starts');

%% M step of DSLDA
model.eta1 = update_eta_fast_NPDSLDA(model.lambda(1:model.K1,:), MaxFun, model.eta1);
disp('eta1 done');
if(option==1)
 model.eta2 = update_eta_fast_NPDSLDA(model.lambda(model.K1+1:end,:), MaxFun, model.eta2);
end
disp('eta2 done');
% value = likelihood_NPDSLDA(model, data, option)
% if (compareval(value, maxvalue))
%     maxvalue = value;
% else
%     error('Incorrect after eta');
% end

model   = update_gammazero_NPDSLDA(model, MaxFun);
disp('gammazero done');
% value = likelihood_NPDSLDA(model, data, option)
% if (compareval(value, maxvalue))
%     maxvalue = value;
% else
%     error('Incorrect after gammazero');
% end

model   = update_alphazero_NPDSLDA(model, MaxFun);
disp('alphazero done');
% value = likelihood_NPDSLDA(model, data, option)
% if (compareval(value, maxvalue))
%     maxvalue = value;
% else
%     error('Incorrect after alphazero');
% end

if(option==1)
    model   = update_uzero_NPDSLDA(model, MaxFun);
    disp('uzero done');
    % value = likelihood_NPDSLDA(model, data, option)
    % if (compareval(value, maxvalue))
    %     maxvalue = value;
    % else
    %     error('Incorrect after uzero');
    % end
end

model  = update_r(model, data, option);

%%   value3 = likelihood_NPDSLDA(model, data, option)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% don't include the following unless you calculate the slack variables from svm package and include that in the lower bound calculation
%    maxvalue = value3;
%     if (compareval(value3, maxvalue))
%         maxvalue = value3;
%     else
%         error('Incorrect after eta');
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('M step done');

end
