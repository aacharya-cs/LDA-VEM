function [trdataova, testdataova] = get_medLDAovadata(trdata, testdata, classind)

N1 = length(trdata.classlabels);
indpositive = find(trdata.classlabels==classind); % positive classes
indnegative = setdiff([1:N1],indpositive);

%% sub-sample the negative class
negind = SelRandomVec(length(indnegative),length(indpositive));
%negind = SelRandomVec(length(indnegative),length(indnegative));
indnegative = indnegative(negind);

for i=1:length(indpositive)
    trdataova.windex{i} = trdata.windex{indpositive(i)};
    trdataova.wcount{i} = trdata.wcount{indpositive(i)};
    trdataova.annotations(i,:) = trdata.annotations(indpositive(i),:);
    trdataova.nwordspdoc(i,1) = trdata.nwordspdoc(indpositive(i));
    trdataova.classlabels(i,1) = 1;
end

for i=1:length(indnegative)
    temp = length(indpositive)+i;
    trdataova.windex{temp} = trdata.windex{indnegative(i)};
    trdataova.wcount{temp} = trdata.wcount{indnegative(i)};
    trdataova.annotations(temp,:) = trdata.annotations(indnegative(i),:);
    trdataova.nwordspdoc(temp,1) = trdata.nwordspdoc(indnegative(i));
    trdataova.classlabels(temp,1) = 2;
end
trdataova.Y  = 2;
trdataova.k1 = trdata.k1;
trdataova.k2 = trdata.k2;
trdataova.V  = trdata.V;

N2 = length(testdata.classlabels);
indpositive = find(testdata.classlabels==classind); % positive classes
indnegative = setdiff([1:N2],indpositive);
testdataova = testdata;
testdataova.classlabels(indpositive) = 1;
testdataova.classlabels(indnegative) = 2;
testdataova.Y = 2;

end
