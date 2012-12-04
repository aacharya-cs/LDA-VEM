function [T] = contractT(Tmat)

Tmatind = find(Tmat'==1);
k = size(Tmat,2);
T = mod(Tmatind,k);
T(find(T==0))=k;

end