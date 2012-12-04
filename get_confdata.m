function [trdata, testdata, selindex, classnames, MCMacc] = get_confdata(pathname, p1, p2, p3, classfilename, k2, minvtopic)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% treats the last file name as the target class
%% p1: percentage or number of instances from the source data in training
%% p2: percentage or number of instances from the target data in training
%% p3: percentage of training data that have missing labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classnames   = textread(classfilename,'%s\n');
for i=1:size(classnames,1)-1
    source_class{i} = classnames{i};
end
target_class = classnames{end};

A = load('vocabList.mat');
trdata.V    = size(A.vocabList,2); %% code book size or vocabulary size
testdata.V  = size(A.vocabList,2);

A = load('tagList.mat');
trdata.k1   = size(A.tagList,2);  %% maximum number of visible topics
testdata.k1 = size(A.tagList,2);

trdata.k2   = k2;  %% number of latent topics
testdata.k2 = k2;

[M N] = size(classnames);

trdata.Y    = max(M,N);  %% number of classes
testdata.Y  = max(M,N);

%% collect training data
count = 1;
for i=1:trdata.Y
    invalidindex = [];
    if(i<trdata.Y)
        filename = [pathname source_class{i} '_train.mat'];
        A = load(filename);
        p = p1;
    else
        filename = [pathname target_class '_train.mat'];
        A = load(filename);
        p = p2;
    end
    N = size(A.wcount,2);
    selindtrain = SelRandomVec(N, round(p*N));
    
    
    %% select only a subset of training data to have labels
    Ntrain = length(selindtrain);
    selindtrlabels = SelRandomVec(Ntrain, round(p3*Ntrain));
    selindtrainbinarylabels = zeros(1,Ntrain);
    selindtrainbinarylabels(selindtrlabels) = 1;
    %% selection done
    
    for j=1:Ntrain
        temp = find(A.annotations(selindtrain(j),:)==1);
        if(length(temp)<minvtopic)
            invalidindex = [invalidindex j];
        else
            trdata.windex{count} = A.windex{selindtrain(j)};
            trdata.wcount{count} = A.wcount{selindtrain(j)};
            trdata.annotations(count,:) = A.annotations(selindtrain(j),:);
            trdata.classlabels(count,1) = i*selindtrainbinarylabels(j); %% assigned zero labels for missing labels
            trdata.nwordspdoc(count,1) = sum(A.wcount{selindtrain(j)});
            count = count + 1;
        end
    end
    selindtrain(invalidindex) = [];
    selindtrlabels(find(selindtrlabels>length(selindtrain))) = [];
    selindex{i}.selindtrain = selindtrain(selindtrlabels);
    %%length(selindex{i}.selindtrain)
end

%% collect test data
maxnum = 0;
count = 1;
for i=1:testdata.Y
    invalidindex = [];
    if(i<testdata.Y)
        filename = [pathname source_class{i} '_test.mat'];
        A = load(filename);
        p = 1;
    else
        filename = [pathname target_class '_test.mat'];
        A = load(filename);
        p = 1;
    end
    
    N = size(A.wcount,2);
    selindtest = SelRandomVec(N, round(p*N));
    for j=1:length(selindtest)
        temp = find(A.annotations(selindtest(j),:)==1);
        if(length(temp)<minvtopic)
            invalidindex = [invalidindex j];
        else
            testdata.windex{count} = A.windex{selindtest(j)};
            testdata.wcount{count} = A.wcount{selindtest(j)};
            testdata.annotations(count,:) = A.annotations(selindtest(j),:);
            testdata.classlabels(count,1) = i;
            testdata.nwordspdoc(count,1) = sum(A.wcount{selindtest(j)});
            count = count + 1;
        end
    end
    selindtest(invalidindex) = [];
    selindex{i}.selindtest  = selindtest;
    
    if(length(selindtest)>maxnum)
        maxnum = length(selindtest);
    end
end

MCMacc = maxnum/length(testdata.classlabels);


%% for comparing with the code of onlineHDP
converttoSLDA(trdata, testdata);

%%getstat(trdata);

end

