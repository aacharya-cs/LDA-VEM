function [testmodel] = inittestmodel(model, data, option, p)

testmodel.N   = size(data.wcount,2);
testmodel.K   = model.K;
testmodel.V   = model.V;

%% model parameters copied from the model obtained in training

if(option>=4)
    testmodel.k1      = data.k1;
    testmodel.k2      = data.k2;
    testmodel.epsilon = model.epsilon;
    testmodel.alpha1  = model.alpha1;
    testmodel.alpha2  = model.alpha2;
else
    testmodel.alpha   = model.alpha;
end
testmodel.log_beta    = model.log_beta;  % dimension K*V

testmodel.MINVALUE = model.MINVALUE;
testmodel.phase    = 0;
testmodel.option   = option;

if(option>=3)
    testmodel.Y   = data.Y;
    testmodel.C2  = model.C2;
    testmodel.C1  = model.C1;
    testmodel.eta = model.eta;   % Y*K -- another model parameter
end

%% variational parameters initialized
testmodel.gamma  = p*ones(testmodel.N,testmodel.K);   % dimension N*K
%% new initialization
tmp1  = ones(testmodel.N, testmodel.K)/testmodel.K;
tmp2  = repmat(data.nwordspdoc,1,testmodel.K)/testmodel.K;
testmodel.gamma    = tmp1 + tmp2;

if(option==7) %%DSLDA-OSST;
    ind = [testmodel.k1+1:testmodel.K];
    testmodel.gamma(:,ind) = 0;
end

if(option==8) %%DSLDA-NSLT2;
    ind = [1:testmodel.k1];
    testmodel.gamma(:,ind) = 0;
end

for n=1:testmodel.N
    temp = ones(max(size(data.windex{n})), testmodel.K); %% uniform initialization of phi's
    temp = temp./repmat(sum(temp,2),1,model.K);
    testmodel.phi{n} = temp;
    if(size(temp,1)==1 && size(temp,2)==1)
        n
        data.windex{n}'
        error('problem in init');
    end
end

end
