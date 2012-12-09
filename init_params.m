function [model] = init_params(data, K, V, p, option, phase)

% Random Initialization of model and variational parameters

randind = 1;

model.N   = size(data.wcount,2);
model.V   = data.V;

if(option>=4)
    model.k1      = data.k1;
    model.k2      = data.k2;
    if(option==7) % DSLDA-OSST
        model.epsilon = 0.999999;
    elseif(option==8) % DSLDA-NSLT2
        model.epsilon = 0.0000000001;
    else
        model.epsilon = 0.80;  %% weight of supervised topics
    end
    model.alpha1  = randind*ones(1,model.k1)/model.k1;   %% the supervised topics are shared across DSLDA and DSLDA with no shared latent topic
    
    if (option==5 || option==8)     %% DSLDA-NSLT-Ray
        model.k2 = round(data.k2/data.Y)*data.Y;
        model.K = model.k1+model.k2;
        model.alpha2  = randind*ones(1,model.k2)/model.k2;
    elseif (option==6) %% DSLDA-NSLT-Ayan
        model.K = model.k1+model.k2*data.Y;
        model.alpha2  = randind*ones(1,model.k2*data.Y);
    else
        model.K = model.k1+model.k2;
        model.alpha2  = randind*ones(1,model.k2)/model.k2;
    end
else
    model.K = data.k1+data.k2;    %% total number of topics are held fixed across models
    model.alpha   = randind*rand(1,model.K)/model.K;
end

model.gamma    = randind*rand(model.N,model.K)*(model.N/model.K);   % dimension N*K   %% have to change

%% new initialization
tmp1  = ones(model.N, model.K)/model.K;
tmp2  = repmat(data.nwordspdoc,1,model.K)/model.K;
model.gamma    = tmp1 + tmp2;

temp1          = randind*rand(model.K,V);
model.log_beta = log(temp1./repmat(sum(temp1,2), 1, model.V));  % dimension K*V
for n=1:model.N
    temp = randind*ones(max(size(data.windex{n})), model.K); %% uniform initialization of phi's
    temp = temp./repmat(sum(temp,2),1,model.K);
    model.phi{n} = temp;
    if(size(temp,1)==1 && size(temp,2)==1)
        n
        data.windex{n}'
        error('problem in init');
    end
end
% cell array of dimension N*1
% nth element is of dimension length(windex{n})*K


model.MINVALUE = 0.0000000000000000000000001;
model.phase    = phase;
model.option   = option;

if(option>=3)
    model.Y   = data.Y;
    model.C2  = 0.1;
    model.C1  = 0.1;
    model.mu  = randind*zeros(model.N,model.Y);   % N*Y  dual variables // zero because we have no idea of the duals right now
    model.eta = randind*rand(model.Y,model.K);   % Y*K  svm weights
end


if(model.option==5 || model.option==8) %% DSLDA-NSLT1, DSLDA-NSLT2
    ind2  = rand(model.N,model.k1);
    ind2  = [ind2 2*rand(model.N,model.K-model.k1)];
    ind21 = repmat(data.classlabels-1,1,(model.k2/data.Y))*(model.k2/data.Y);
    ind22 = repmat([1:(model.k2/data.Y)],model.N,1);
    ind23 = ind21+ind22+model.k1;
    ind24 = (ind23-1)*model.N + repmat([1:model.N]',1,size(ind23,2));
    ind2(ind24(:)) = 0;
    ind3 = find(ind2==2);
    model.gamma(ind3) = 0; % zero out indices which are not active among latent topics
end

if(model.option==6) %% DSLDA-NSLT-Ayan
    ind2  = rand(model.N,model.k1);
    ind2  = [ind2 2*rand(model.N,model.K-model.k1)];
    ind21 = repmat(data.classlabels-1,1, model.k2)*model.k2;
    ind22 = repmat([1:model.k2],model.N,1);
    ind23 = ind21+ind22+model.k1;
    ind24 = (ind23-1)*model.N + repmat([1:model.N]',1,size(ind23,2));
    ind2(ind24(:)) = 0;
    ind3 = find(ind2==2);
    model.gamma(ind3) = 0; % zero out indices which are not active among latent topics
end

%% initialization of sufficient statistics
%% random
model.ss_topicword = 10*ones(model.K, model.V) + rand([model.K, model.V]);
model.ss_topic = sum(model.ss_topicword,2)';
%% initialization of beta based on sufficient statistics
model.log_beta = update_beta_cpp (model);

end
