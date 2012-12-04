function [vec3] = lower_partial_sum(vec1)

%% computes \sum_{n'=1}^{n-1}vec1_{n'}{:}
%% vec1 should be a column vector or a matrix 

vec2 = ones(1,size(vec1,1));
vec2(1) = 0;
vec2 = tril(toeplitz(vec2));
vec3 = vec2*vec1;

end