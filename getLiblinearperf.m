function [svmperf, trdata, trlabel, testdata, testlabel] = getLiblinearperf(pathname, selindex, classnames, optionval, cval)

trdata = [];
testdata = [];
trlabel = [];
testlabel = [];
[M N] = size(classnames);

classnames

for i=1:max(M,N)
    trfile = [pathname classnames{i} '_train.mat'];
    A = load(trfile);
    testfile = [pathname classnames{i} '_test.mat'];
    B = load(testfile);
    
    trdata    = [trdata; A.liblindata(selindex{i}.selindtrain,:)];
    trlabel   = [trlabel; i*ones(length(selindex{i}.selindtrain),1)];
    testdata  = [testdata; B.liblindata(selindex{i}.selindtest,:)];
    testlabel = [testlabel; i*ones(length(selindex{i}.selindtest),1)];
end


size(trdata)

option = [' -s ' num2str(optionval) ' -c ' num2str(cval) ' -e 0.1'];
model = train(trlabel, trdata, option);

[predclassestrain, ~, ~] = predict(trlabel, trdata, model);
[predclassestest, ~, ~]  = predict(testlabel, testdata, model);

[F_score_train, confmat_train] = cal_F_score(trlabel, predclassestrain);
[F_score_test, confmat_test]   = cal_F_score(testlabel, predclassestest);

svmperf.F_score_train = F_score_train;
svmperf.F_score_test  = F_score_test;
svmperf.confmat_train = confmat_train;
svmperf.confmat_test  = confmat_test;
svmperf.binacc = binary_accuracy(confmat_test);

%%pause

end
