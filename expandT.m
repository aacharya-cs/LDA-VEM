function [Tex] = expandT(T,Clnum)

testsize=size(T,1);
Tex = zeros(Clnum,testsize);
temp = T+Clnum*[0:testsize-1]';
Tex(temp)=1;
Tex = Tex';

end
