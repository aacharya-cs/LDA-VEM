function [Elogsticks] = E_logsticks(sticks)

dig_sum = psi(sum(sticks, 1));
ElogW   = psi(sticks(1,:)) - dig_sum;
Elog1_W = psi(sticks(2,:)) - dig_sum;

n = size(sticks,2) + 1;
Elogsticks = zeros(1,n);
Elogsticks(1:n-1) = ElogW;
Elogsticks(2:end) = Elogsticks(2:end) + cumsum(Elog1_W);


end
