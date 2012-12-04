function [F_score, confmat] = cal_F_score(A, B)

%% returns F-score and the the confusion matrix

confmat = confusionmat(A,B);

if(length(unique(A))>2)       %% for multi-class problems  
    rowsum  = sum(confmat,2);
    colsum  = (sum(confmat,1))';
    precisionVal = diag(confmat)./rowsum;  %% either one of these would be precision or recall -- don't care which is what now
    recallVal    = diag(confmat)./colsum;   
elseif(length(unique(A))==2)  %% for binary problems  
    rowsum  = sum(confmat(1,:),2);
    colsum  = (sum(confmat(:,1),1))';  
    precisionVal = confmat(1,1)/rowsum;  %% either one of these would be precision or recall -- don't care which is what now
    recallVal    = confmat(1,1)/colsum;
else
    error('only one class in the groundtruth');
end

F_score = 2*(precisionVal.*recallVal)./(precisionVal+recallVal);



end
