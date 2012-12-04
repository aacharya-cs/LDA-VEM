function [model] = M_step(model, data, MaxFun, maxvalue)

disp('M step starts');

%% M step of DSLDA
model.eta = update_eta_fast_NPDSLDA(model.lambda, MaxFun, model.eta);
disp('eta done');
% value = likelihood_NPDSLDA(model, data)
% if (compareval(value, maxvalue))
%     maxvalue = value;
% else
%     error('Incorrect after eta');
% end

model   = update_alphazero_NPDSLDA(model, MaxFun);
disp('alphazero done');
% value = likelihood_NPDSLDA(model, data)
% if (compareval(value, maxvalue))
%     maxvalue = value;
% else
%     error('Incorrect after alphazero');
% end

model   = update_uzero_NPDSLDA(model, MaxFun);
disp('uzero done');
% value = likelihood_NPDSLDA(model, data)
% if (compareval(value, maxvalue))
%     maxvalue = value;
% else
%     error('Incorrect after uzero');
% end

model   = update_gammazero_NPDSLDA(model, MaxFun);
disp('gammazero done');
% value = likelihood_NPDSLDA(model, data)
% if (compareval(value, maxvalue))
%     maxvalue = value;
% else
%     error('Incorrect after gammazero');
% end

model  = update_r(model, data, model.sumzeta);

%%   value3 = likelihood_NPDSLDA(model, data)

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
