function [vec3] = partial_sum(A)

%% if A is the input matrix, then the output is
%% A(:,t) = \sum_{k=t}^{T} A(:,k)
%% any 1D vector should be a row vector; the returned vector is also a row vector

vec1 = flipdim(A,2);
vec2 = cumsum(vec1,2);
vec3 = flipdim(vec2,2);

end