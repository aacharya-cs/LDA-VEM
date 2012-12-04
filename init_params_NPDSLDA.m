function [model, data] = init_params_NPDSLDA(data, K1, T)

%% random Initialization of model and variational parameters
model.N        = size(data.wcount,2);
model.V        = data.V;
model.T        = T;
model.K1       = K1;
model.K2       = data.k1; %% number of supervised and latent topics flipped
model.epsilon  = 0.8;    %% weight of supervised topics
model.phase    = 1;
model.C1       = 1;
model.C2       = 5;
model.MINVALUE = 1e-25;

%% model parameters
model.gammazero = [1 2];  %% prior over global sticks
model.alphazero = [1 2];  %% prior over local sticks
model.uzero     = rand(1,model.K2);
model.eta       = ones(1,model.V);    %% dimension 1*V
model.r         = ones(data.Y,(model.T+model.K2));              %% dimension Y*(T+K2)
model.dmu       = zeros(model.N,data.Y); %% N*Y  dual variables // zero because we have no idea of the duals right now

%% variational parameters
model.u         = model.gammazero(1)*ones(1,model.K1-1);        %% dimension 1*(K1-1)
model.v         = model.gammazero(2)*ones(1,model.K1-1);        %% dimension 1*(K1-1)
model.a         = model.alphazero(1)*ones(model.N,model.T-1);   %% dimension N*(T-1)
model.b         = model.alphazero(2)*ones(model.N,model.T-1);   %% dimension N*(T-1)
model.lambda    = rand(model.K1+model.K2,model.V);   % dimension (K1+K2)*V
model.mun       = rand(model.N,model.K2);  % dimension N*K2
model.smallphi  = rand(model.N, model.T, model.K1); % dimension N*T*K1
if(model.phase==1)
    model.mun = model.mun.*data.annotations;
end

%% normalization of variational parameters
for n=1:model.N
    %% smallphi
    temp1 = squeeze(model.smallphi(n,:,:));
    temp2 = sum(temp1,2);
    temp1 = temp1./repmat(temp2,1,size(temp1,2));
    model.smallphi(n,:,:) = temp1;
    
    %% zeta
    temp3 = rand(max(size(data.windex{n})), model.T);  
    temp4 = rand(max(size(data.windex{n})), model.K2); 
    if(model.phase==1)
        temp4 = temp4.*repmat(data.annotations(n,:),max(size(data.windex{n})),1);
    end
    temp5  = [temp3 temp4];
    temp5  = temp5./repmat(sum(temp5,2),1,(model.T+model.K2));
    model.zeta{n} = temp5;
    %% for checking errors
    if(size(temp5,1)==1 && size(temp5,2)==1)
        n
        data.windex{n}'
        error('problem in init');
    end
end
model.sumzeta  = sum_zeta(model, data);

end

