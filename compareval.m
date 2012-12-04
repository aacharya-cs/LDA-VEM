function [boolnum] = compareval(a, b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

temp = (abs(a-b)/abs(a));
if(a>=b ||temp<0.00000001)
    boolnum=1;
else
    boolnum=0;
end

end

