function [] = getstat(trdata)

[N vtopics] = size(trdata.annotations);
classnum    = length(unique(trdata.classlabels));
totwords  = [];
tottopics = [];
classlabels = [];
topicdoccount = zeros(1,vtopics);
H = zeros(classnum,vtopics);

for i=1:N
    %  totwords  = [totwords trdata.nwordspdoc(i)];
    %  tottopics = [tottopics length(find(trdata.annotations(i,:)==1))];
    H(trdata.classlabels(i),:) = H(trdata.classlabels(i),:) + trdata.annotations(i,:);
    %  topicdoccount(1,ind) = topicdoccount(1,ind) + 1;
    classlabels = [classlabels trdata.classlabels(i)];
end

uniqclass = unique(classlabels)
C = sort(uniqclass,'ascend');
B = hist(classlabels, C);
H = H./repmat(B',1,vtopics);

A = zeros(classnum,classnum);
for i=1:classnum
    for j=i+1:classnum
        temp = sum([H(i,:)-H(j,:)].^2);
        A(i,j) = sqrt(temp);
        A(j,i) = A(i,j);
    end
end

H
A = A/max(max(A));
S = (1-A);
% %hist(totwords);
% %xlabel('number of words', 'FontSize', 14);
% %ylabel('number of documents', 'FontSize', 14);
% %title('distribution of words per document', 'FontSize', 14);
% %h = findobj(gca,'Type','patch');
% %set(h,'FaceColor','r','EdgeColor','w');
%
% mean(totwords);
% std(totwords);
%
% %hist(tottopics);
% %xlabel('number of visible topics', 'FontSize', 14);
% %ylabel('number of documents', 'FontSize', 14);
% %title('distribution of visible topics per document', 'FontSize', 14);
% %h = findobj(gca,'Type','patch');
% %set(h,'FaceColor','r','EdgeColor','w');
%
%
% mean(tottopics);
% std(tottopics);
%
% topicdoccount;
%
% %bar(sort(topicdoccount,'descend'))
% %xlabel('indices of visible topics', 'FontSize', 14);
% %ylabel('number of documents', 'FontSize', 14);
% %title('distribution of number of documents per visible topic', 'FontSize', 14);
% %h = findobj(gca,'Type','patch');
% %set(h,'FaceColor','r','EdgeColor','w');
%
% ind = hist(classlabels);
% ind(find(ind)==0) = [];
% ind = sort(ind,'descend');
% bar(ind(1:classnum));
%
% xlabel('indices of classes', 'FontSize', 14);
% ylabel('number of documents', 'FontSize', 14);
% title('distribution of number of documents per class', 'FontSize', 14);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w');

pause

end
