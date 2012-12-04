function [trdata] = createsynthtrdata()

N = 5;
for i=1:N
  temp = rand(1,4);
  temp(find(temp>0.5))  = 1;
  temp(find(temp<=0.5)) = 0;
  trdata.annotations(i,:)= temp;
  trdata.classlabels(i) = ceil(2*rand);
end

end
