function [testmodel] = inittestmodel_NPDSLDA(trmodel, data)

testmodel.N         = size(data.wcount,2);
testmodel.V         = data.V;
testmodel.T         = trmodel.T;
testmodel.K1        = trmodel.K1;
testmodel.K2        = trmodel.K2;         %% number of supervised and latent topics flipped
testmodel.epsilon   = trmodel.epsilon;    %% weight of supervised topics
testmodel.phase     = 0;
testmodel.C1        = trmodel.C1;
testmodel.C2        = trmodel.C2;
testmodel.MINVALUE  = 1e-25;

%% model parameters
testmodel.gammazero = trmodel.gammazero;
testmodel.alphazero = trmodel.alphazero;
testmodel.uzero     = trmodel.uzero;
testmodel.eta       = trmodel.eta;    %% dimension 1*V
testmodel.r         = trmodel.r;      %% dimension Y*(T+K2)

%% global variational parameters
testmodel.u         = trmodel.u;        %% dimension 1*(K1-1)
testmodel.v         = trmodel.v;        %% dimension 1*(K1-1)
testmodel.lambda    = trmodel.lambda;   %% dimension (K1+K2)*V

%% local variational parameters
testmodel.a         = rand(testmodel.N,testmodel.T);   %% dimension N*T
testmodel.b         = rand(testmodel.N,testmodel.T);   %% dimension N*T
testmodel.mun       = rand(testmodel.N,testmodel.K2);  %% dimension N*K2
testmodel.smallphi  = rand(testmodel.N, testmodel.T, testmodel.K1); %% dimension N*T*K1

%% normalization of local variational parameters
for n=1:testmodel.N
    temp1 = squeeze(testmodel.smallphi(n,:,:));
    temp2 = sum(temp1,2);
    temp1 = temp1./repmat(temp2,1,size(temp1,2));
    testmodel.smallphi(n,:,:) = temp1;
    
    temp3 = rand(max(size(data.windex{n})), testmodel.T);  
    temp4 = rand(max(size(data.windex{n})), testmodel.K2); 
    temp5 = [temp3 temp4];
    temp5 = temp5./repmat(sum(temp5,2),1,(testmodel.T+testmodel.K2));
    testmodel.zeta{n} = temp5;
    %% for checking errors
    if(size(temp5,1)==1 && size(temp5,2)==1)
        n
        data.windex{n}'
        error('problem in init');
    end   
end
testmodel.sumzeta  = sum_zeta(testmodel, data);

end
