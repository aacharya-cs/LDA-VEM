function [F_score, confmat, probassign, binacc] = cal_accuracy_NPDSLDA(model, data, option)

probassign = [];
temp  = model.ss_features;
if(option==1)
    temp  = temp./repmat(data.nwordspdoc,1,(model.K1+model.K2));
else
    temp  = temp./repmat(data.nwordspdoc,1,model.K1);
end
temp2 = temp*model.r'; %% N*Y
[~,ind] = max(temp2,[],2);

temp3 = exp(temp2)./repmat(sum(exp(temp2),2),1,data.Y);
probassign = temp3(:,1);
[F_score, confmat] = cal_F_score(data.classlabels, ind);
binacc = binary_accuracy(confmat);

end
