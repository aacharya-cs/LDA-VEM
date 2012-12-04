function [ind]  = SelRandomVec(N,k)

%% N: size of vector, k: number of elements to be chosen
%% if k is zero, return the first element

if(k==N)
    ind = [1:N];
elseif (k==0)
    ind = 1;    
else
    permvect = randperm(N);
    ind      = permvect(1:k);
    ind      = sort(ind); 
end

end