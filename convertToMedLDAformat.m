function [strtext] = convertToMedLDAformat(wind, wcount, class)

%% A: wordindex
%% B: wordcounts
%% C: class label

strtext = [num2str(length(wind)) ' ' num2str(class) ' '];
A = [wind' wcount'];
A = A';
A = A(:)';
tempstr1 = num2str(A);
tempstr2 = regexprep(tempstr1,'\s*', ' ');
spaceind = regexp(tempstr2,' ');
oddind   = find(rem([1:length(spaceind)],2)>0);
tempstr2(spaceind(oddind)) = ':';
strtext = [strtext tempstr2];

end
