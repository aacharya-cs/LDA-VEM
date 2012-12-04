function [model] = update_gamma(model, data)

if(model.option>=4)
    alphaterm1 = [repmat(model.alpha1, model.N, 1) zeros(model.N, model.K-model.k1)];
    alphaterm2 = [zeros(model.N, model.k1) repmat(model.alpha2, model.N, 1)];
    term1      = alphaterm1 + alphaterm2;
else
    term1 = repmat(model.alpha, model.N, 1);
end

term2 = sum_phi(model, data);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%% inefficient matlab code..use when the C++ code is not correct
% % term22 = zeros(model.N,model.K);
% %  for n=1:model.N
% %   term2temp = repmat(data.wcount{n}',1,model.K).*model.phi{n};
% %   term22(n,:) = term22(n,:) + sum(term2temp,1);
% %  end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % term2-term22
% % ind1 = find(term1<0);
% % ind2 = find(term2<0);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.gamma = term1+term2;

if (model.phase==1 && (model.option==2||model.option>=4)) % only in training phase and for option labeled LDA and my model
    ind1 = find(data.annotations==0); % since latent topics are appended as extra columns; does not make a difference in indexing
    model.gamma(ind1) = 0; % zero out indices which are not active among supervised topics
    if(model.option==5 || model.option==8) %% DSLDA-NSLT1,DSLDA-NSLT2 
        ind2  = ones(model.N,model.k1);
        ind2  = [ind2 2*ones(model.N,model.K-model.k1)];
        ind21 = repmat(data.classlabels-1,1,(model.k2/data.Y))*(model.k2/data.Y);
        ind22 = repmat([1:(model.k2/data.Y)],model.N,1);
        ind23 = ind21+ind22+model.k1;
        ind24 = (ind23-1)*model.N + repmat([1:model.N]',1,size(ind23,2));
        ind2(ind24(:)) = 0;
        ind3 = find(ind2==2);
        model.gamma(ind3) = 0; % zero out indices which are not active among latent topics
    end
    if(model.option==6) %% DSLDA-NSLT-Ayan
        ind2  = ones(model.N,model.k1);
        ind2  = [ind2 2*ones(model.N,model.K-model.k1)];
        ind21 = repmat(data.classlabels-1,1, model.k2)*model.k2;
        ind22 = repmat([1:model.k2],model.N,1);
        ind23 = ind21+ind22+model.k1;
        ind24 = (ind23-1)*model.N + repmat([1:model.N]',1,size(ind23,2));
        ind2(ind24(:)) = 0;
        ind3 = find(ind2==2);
        model.gamma(ind3) = 0; % zero out indices which are not active among latent topics
    end
end

% % if(model.option==7)
% %     %% DSLDA-OSST; have to do it both in training and test phase; otherwise
% %     %% non-zero gammas (which are supposed to be zero) affect inference in test phase
% %     ind3 = [model.k1+1:model.k2];
% %     model.gamma(:,ind3) = 0; %% zero out indices corresponding to latent topics
% % end
% % 
% % if(model.option==8)
% %     %% DSLDA-OSST; have to do it both in training and test phase; otherwise
% %     %% non-zero gammas (which are supposed to be zero) affect inference in test phase
% %     ind3 = [1:model.k1];
% %     model.gamma(:,ind3) = 0; %% zero out indices corresponding to supervised topics
% % end


end
