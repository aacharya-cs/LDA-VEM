function [model] = update_r(model, data, sumzeta)

cval         = model.C2/model.C1;
option       = [' -s 4'  ' -c ' num2str(cval) ' -e 0.1'];
trdata       = sumzeta./repmat(data.nwordspdoc,1,(model.T+model.K2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % if(model.option~=5)
% %     trdata       = temp./repmat(data.nwordspdoc,1,model.K);
% % else
% %     tempmat  = expandT(data.classlabels,model.Y);
% %     tempmat  = repmat(tempmat,1,model.Y);
% %     ind      = reshape([1:1:(model.Y)^2],model.Y,model.Y);
% %     ind      = ind';
% %     ind      = (ind(:))';
% %     tempmat  = tempmat(:,ind);
% %     
% %     temp1    = temp(:,1:model.k1);  %% shared weights
% %     temp2    = temp(:,model.k1+1:model.K2).*tempmat; %% weights for non-shared latent topics..zero padding is used for "other" classes
% %     temp     = [temp1 temp2];
% %     trdata   = temp./repmat(data.nwordspdoc,1,model.k1+model.k2*model.Y);
% % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% for training only with labeled data
ind = find(data.classlabels>0);
svmtrlabels = data.classlabels(ind);
svmtrdata = sparse(trdata(ind,:));
%% selection done

trainedmodel = train(svmtrlabels, svmtrdata, option);
[~, acc, ~]  = predict(svmtrlabels, svmtrdata, trainedmodel); % test the training data
acc
model.r      = trainedmodel.w;
model.dmu    = trainedmodel.alpha;

end
