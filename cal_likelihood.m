function [value] = cal_likelihood (model, data)

%%%% modification of terms for labeled LDA or my model

ind  = [];
if(model.phase==1 && (model.option==2||model.option>=4))
    if(model.option==2)
        
        alphaterm = repmat(model.alpha, model.N, 1);
        M = [data.annotations];
        ind = find(M==0);
        mod_alpha  = (M.*alphaterm);   %alpha_{ni}  -- we need to ignore co-ordinates which are zero
        
        gammaln_mod_alpha = gammaln(mod_alpha);
        gammaln_mod_alpha(find(gammaln_mod_alpha==inf)) = 0;               %gammaln(alpha_{ni}) inf replaced by zero as we need to sum over the co-ordinates
        
        gammaln_gamma  = gammaln(model.gamma);
        gammaln_gamma(find(gammaln_gamma==inf)) = 0;      % some gammas are zeroed out in training phase, need to undo them here
        
        psi_gamma  = psi(model.gamma);
        psi_gamma(find(psi_gamma==-inf)) = 0;      % some gammas are zeroed out in training phase, need to undo them here
        
    end
    if(model.option>=4)
        
        M = [data.annotations ones(model.N, model.K-model.k1)];
        ind = find(M==0);
        alphaterm1 = [repmat(model.alpha1, model.N, 1) zeros(model.N, model.K-model.k1)];
        alphaterm2 = [zeros(model.N, model.k1) repmat(model.alpha2, model.N, 1)];
        
        if(model.option==5)
            ind2  = ones(model.N,model.k1);
            ind2  = [ind2 2*ones(model.N,model.K-model.k1)];
            ind21 = repmat(data.classlabels-1,1,(model.k2/data.Y))*(model.k2/data.Y);
            ind22 = repmat([1:(model.k2/data.Y)],model.N,1);
            ind23 = ind21+ind22+model.k1;
            ind24 = (ind23-1)*model.N + repmat([1:model.N]',1,size(ind23,2));
            ind2(ind24(:)) = 0;
            ind3 = find(ind2==2);
            alphaterm2(ind3) = 0; % zero out indices which are not active among latent topics
            ind = [ind; ind3];
        end
        if(model.option==6)
            ind2  = ones(model.N,model.k1);
            ind2  = [ind2 2*ones(model.N,model.K-model.k1)];
            ind21 = repmat(data.classlabels-1,1,model.k2)*model.k2;
            ind22 = repmat([1:model.k2],model.N,1);
            ind23 = ind21+ind22+model.k1;
            ind24 = (ind23-1)*model.N + repmat([1:model.N]',1,size(ind23,2));
            ind2(ind24(:)) = 0;
            ind3 = find(ind2==2);
            alphaterm2(ind3) = 0; % zero out indices which are not active among latent topics
            ind = [ind; ind3];
        end
        
        
        mod_alpha  = (M.*alphaterm1)*model.epsilon + (1-model.epsilon)*alphaterm2;   %alpha_{ni}  -- we need to ignore co-ordinates which are zero
        gammaln_mod_alpha = gammaln(mod_alpha);
        gammaln_mod_alpha(find(gammaln_mod_alpha==inf)) = 0;               %gammaln(alpha_{ni}) inf replaced by zero as we need to sum over the co-ordinates
        
        gammaln_gamma  = gammaln(model.gamma);
        gammaln_gamma(find(gammaln_gamma==inf)) = 0;      % some gammas are zeroed out in training phase, need to undo them here
        
        psi_gamma  = psi(model.gamma);
        psi_gamma(find(psi_gamma==-inf)) = 0;      % some gammas are zeroed out in training phase, need to undo them here
        
    end
else
    if(model.option>=4)
        alphaterm1 = [repmat(model.alpha1, model.N, 1) zeros(model.N, model.K-model.k1)];
        alphaterm2 = [zeros(model.N, model.k1) repmat(model.alpha2, model.N, 1)];
        mod_alpha  = alphaterm1 + alphaterm2;
    else
        mod_alpha  = repmat(model.alpha,model.N,1);
    end
    
    gammaln_mod_alpha = gammaln(mod_alpha);
    gammaln_mod_alpha(find(gammaln_mod_alpha==inf)) = 0;
    gammaln_gamma  = gammaln(model.gamma);
    gammaln_gamma(find(gammaln_gamma==inf)) = 0;
    psi_gamma  = psi(model.gamma);
    psi_gamma(find(psi_gamma==-inf)) = 0;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% term1
term11 = sum(gammaln(sum(mod_alpha,2))) -sum(sum(gammaln_mod_alpha));
term12 = (mod_alpha-1).*[psi_gamma - repmat(psi(sum(model.gamma,2)),1,model.K)];
term12(ind) = 0;
term12 = sum(sum(term12));
term1  = (term11+term12);

%% term2
term21 = sum_phi(model, data);
term22 = [psi_gamma-repmat(psi(sum(model.gamma,2)),1,model.K)];
term23 = sum(sum(term21.*term22));
term24 = 0;
if(model.option>=4)
    temp1 = term21(:,1:model.k1)*log(model.epsilon);
    temp2 = term21(:,model.k1+1:model.K)*log(1-model.epsilon);
    temp3 = [temp1 temp2];
    term24 = sum(sum(temp3));
end
term2  = term23 + term24;

%% term3
term3 = sum(sum((model.ss_topicword).*model.log_beta));

%% term4
term41 = sum([gammaln(sum(model.gamma,2))-sum(gammaln_gamma,2)]);
term42 = (model.gamma-1).*[psi_gamma - repmat(psi(sum(model.gamma,2)),1,model.K)];
term42(ind) = 0;
term4  = -(term41+sum(sum(term42)));

%% term5
term5 = 0;
for n=1:model.N
    term51 = repmat(data.wcount{n}',1,model.K).*model.phi{n}.*log(model.phi{n});
    ind2   = find(isnan(term51)==1);
    term51(ind2) = 0;
    term5  = term5 + sum(sum(term51));
end
term5 = -term5;

term6 = 0;
if(model.option>=3)
    term61 = sum_phi(model, data);
    for n=1:model.N
        term62 = repmat(model.eta(data.classlabels(n),:),model.Y,1);
        if((size(term62,1)~=size(model.eta,1)) ||(size(term62,2)~=size(model.eta,2)))
            error('problem');
        end
        if(data.nwordspdoc(n)==0)
            n
            error('no word problem');
        end
        if(model.phase==1)
            term6  = term6 + [1/data.nwordspdoc(n)]*sum([model.mu(n,:)*(term62-model.eta)].*term61(n,:));
        end
    end
end

%% final value
value = term1+term2+term3+term4+term5+term6;

if(isnan(value) || value==inf || value==-inf)
    keyboard;
end

%disp('inside cal_likelihood');

end
