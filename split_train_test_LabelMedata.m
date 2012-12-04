function [] = split_train_test_LabelMedata(sourcepath, savepath)


confname = {'Coast', 'Forest', 'Highway', 'Insidecity', 'Mountain', 'Opencountry', 'Street', 'Tallbuilding'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:size(confname,2)
    
    confname{ii}
    datafilename = [sourcepath confname{ii} '_LabelMe.mat'];    
    savetrdata = [savepath confname{ii} '_train.mat'];
    savetestdata = [savepath confname{ii} '_test.mat'];
    
    A = load(datafilename);
    N = size(A.vwordscount,2);  % total number of examples in the current class
    
    trind  = SelRandomVec(N,round(N/2));
    testind = setdiff([1:N],trind);
    
    wcount = [];
    windex = [];
    nwordspdoc = [];
    class_labels = [];
    annotations = [];
    liblindata = [];
    
    for i=1:length(trind)
        wcount{i} = A.vwordscount{trind(i)};
        windex{i} = A.vwordsindex{trind(i)};
        nwordspdoc(i) = sum(A.vwordscount{trind(i)});
        class_labels(i) = ii;
        annotations(i,:) = A.annotations(trind(i),:);
    end
    liblindata = A.liblindata(trind,:);
    save(savetrdata, 'wcount', 'windex', 'nwordspdoc', 'class_labels', 'annotations', 'liblindata');
    
    wcount = [];
    windex = [];
    nwordspdoc = [];
    class_labels = [];
    annotations = [];
    liblindata = [];
    
    for i=1:length(testind)
        wcount{i} = A.vwordscount{testind(i)};
        windex{i} = A.vwordsindex{testind(i)};
        nwordspdoc(i) = sum(A.vwordscount{testind(i)});
        class_labels(i) = ii;
        annotations(i,:) = A.annotations(testind(i),:);
    end
    liblindata = A.liblindata(testind,:);
    save(savetestdata,'wcount', 'windex', 'nwordspdoc', 'class_labels', 'annotations', 'liblindata');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
