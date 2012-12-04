function [F_score, confmat, probassign, binacc] = cal_accuracy(model, data)

probassign = [];
temp  = sum_phi(model, data);
temp  = temp./repmat(data.nwordspdoc,1,model.K);
temp2 = temp*model.eta'; %% N*Y
[~,ind] = max(temp2,[],2);

if(model.option==3)
 temp3 = exp(temp2)./repmat(sum(exp(temp2),2),1,model.Y);
 probassign = temp3(:,1);
end

[F_score, confmat] = cal_F_score(data.classlabels, ind);
binacc = binary_accuracy(confmat);

end
