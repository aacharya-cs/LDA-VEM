function [strtext] = convertTosLDAformat(wind, wcount)

%% A: wordindex
%% B: wordcounts

wind = wind - 1; %% indexing in sLDA starts from zero
strtext = [num2str(length(wind)) ' '];
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
