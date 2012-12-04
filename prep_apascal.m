function [trdata, testdata, selindex, classnames] = prep_apascal(pathname, p1, p2, classfilename, k2, minvtopic)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% treats the last file name as the target class
%% p1: percentage or number of instances from the source data in training
%% p2: percentage or number of instances from the target data in training
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathname1 = pathname;
pathname2 = pathname;

classnames   = textread(classfilename,'%s\n');
for i=1:size(classnames,1)-1
    source_class{i} = classnames{i};
end
target_class = classnames{end};

trdata.k2   = k2;  %% number of latent topics
testdata.k2 = k2;

[M N] = size(classnames);
trdata.Y    = max(M,N);  %% number of classes
testdata.Y  = max(M,N);

trdata.V    = 500; %% code book size or vocabulary size
testdata.V  = 500;

trdata.k1   = 64;  %% maximum number of visible topics
testdata.k1 = 64;

%% collect training data
count = 1;
for i=1:trdata.Y
    invalidindex = [];
    if(i<trdata.Y)
        filename = [pathname1 source_class{i} '_train.mat'];
        A = load(filename);
        p = p1;
    else
        filename = [pathname1 target_class '_train.mat'];
        A = load(filename);
        p = p2;
    end
    N = size(A.wcount,2);
    selindtrain = SelRandomVec(N, round(p*N));
    for j=1:length(selindtrain)
        if(sum(A.annotations(selindtrain(j),:))==0)
            invalidindex = [invalidindex j];
        else
            if(sum(A.wcount{selindtrain(j)})==0)
            else
                trdata.windex{count} = A.windex{selindtrain(j)};
                trdata.wcount{count} = A.wcount{selindtrain(j)};
                trdata.annotations(count,:) = A.annotations(selindtrain(j),:);
                trdata.classlabels(count,1) = i;
                trdata.nwordspdoc(count,1) = sum(A.wcount{selindtrain(j)});
                count = count + 1; %% to avoid docs without words
                %%error('problem in words');
            end
        end
    end
    if(~isempty(invalidindex))
        %% remove images with no visible topic
        selindtrain(invalidindex) = [];
    end
    selindex{i}.selindtrain = selindtrain;
end

%% collect test data
count = 1;
for i=1:testdata.Y
    invalidindex = [];
    if(i<testdata.Y)
        filename = [pathname2 source_class{i} '_test.mat'];
        A = load(filename);
        p = 1;
    else
        filename = [pathname2 target_class '_test.mat'];
        A = load(filename);
        p = 1;
    end
    
    N = size(A.wcount,2);
    selindtest = SelRandomVec(N, round(p*N));
    for j=1:length(selindtest)
        if(sum(A.annotations(selindtest(j),:))==0)
            invalidindex = [invalidindex j];
        else
            if(sum(A.wcount{selindtest(j)})==0)
            else
                testdata.windex{count} = A.windex{selindtest(j)};
                testdata.wcount{count} = A.wcount{selindtest(j)};
                testdata.annotations(count,:) = A.annotations(selindtest(j),:);
                testdata.classlabels(count,1) = i;
                testdata.nwordspdoc(count,1) = sum(A.wcount{selindtest(j)});
                count = count + 1; %% to avoid docs without words
                %%error('problem in words');
            end
        end
    end
    if(~isempty(invalidindex))
        %% remove images with no visible topic
        selindtest(invalidindex) = [];
    end
    selindex{i}.selindtest  = selindtest;
end

%% getstat(trdata)

end

