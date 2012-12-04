function [tempstr] =  convertToLibSVMformat(matA, label)

temp = full(matA);
ind  = find(temp>0);
A = [ind' (temp(ind))'];
A = A';
B = A(:)'
tempstr1 = num2str(B);
tempstr2 = regexprep(tempstr1,'\s*', ' ');
spaceind = regexp(tempstr2,' ');
oddind   = find(rem([1:length(spaceind)],2)>0);
tempstr2(spaceind(oddind)) = ':'
tempstr = [num2str(length(ind)) ' ' num2str(label) ' ' tempstr2];

end
