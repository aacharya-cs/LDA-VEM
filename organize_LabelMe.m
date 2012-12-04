function [] =  organize_LabelMe()

%% separates the visual words and annotations into categories
classnames = {'Coast', 'Forest', 'Highway', 'Insidecity', 'Mountain', 'Opencountry', 'Street', 'Tallbuilding'};

N = max(size(classnames));
V = 200;

GlobalAttrList = load('GlobalAttrList.mat');
maxvtopics = size(GlobalAttrList.GlobalAttrList,2);

for i=1:N

 filename1 = [classnames{i} '_LabelMe_annotations.mat'];
 filename2 = [classnames{i} '_LabelMe_Vwords.mat'];
 filename3 = [classnames{i} '_LabelMe.mat'];

 A = load(filename1);
 B = load(filename2);
 B = B.A;

 filename = A.FileName;
 M = max(size(filename));
 annotations = zeros(M,maxvtopics);

 vwordsindex = B.vwordsindex; 
 vwordscount = B.vwordscount;  
 liblindata = [];
 for j=1:M
  tempmat = zeros(1,V);
  tempmat(vwordsindex{j}) = vwordscount{j};
  liblindata = [liblindata; sparse(tempmat)];
  annotations(j, A.AttrIndex{j}) = 1;
 end

 save(filename3, 'filename', 'annotations', 'vwordsindex', 'vwordscount', 'liblindata');

end

end
