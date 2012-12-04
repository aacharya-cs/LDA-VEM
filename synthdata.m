function [data] = synthdata(N, V, K, Y, option)

w = rand(N,V);
ind  = find(w<0.2);
w(ind) = 0;
ind  = find(w>=0.2);
w(ind) = 1;
data.nwordspdoc = sum(w,2);

N = size(w,1);
for n=1:N
    ind = find(w(n,:)==1);
    ind = sort(ind);  
    data.windex{n} = ind;
    data.wcount{n} = round(100*rand(1,length(ind)))+1;
end
data.V = size(w,2);

if(option==4)
    k1 = ceil(K/2);
else
    k1 = K;
end

data.annotations = rand(N,k1);
ind  = find(data.annotations<0.2);
data.annotations(ind) = 0;
ind  = find(data.annotations>=0.2);
data.annotations(ind) = 1;

temp = ceil(Y*rand(1,N));
data.classlabels = sort(temp);
data.classlabels = (data.classlabels)'; 

for y = 1:Y
    if(isempty(find(temp==y)))
        error('not all Ys present in training data..try anotther simulation');
    end
end    

data.Y  = Y;
data.k1 = k1;
data.k2 = (K-data.k1);

data.ss_topicword = [];
data.ss_topic     = [];

end
