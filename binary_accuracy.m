function [acctab] = binary_accuracy(A)
% implements binary accuracy from a multi-class confusion matrix
% @ Ayan Acharya, Date 22 Feb, 2012

N      = sum(sum(A));
colsum = sum(A,1)';
rowsum = sum(A,2);
acctab = [N - rowsum - colsum + 2*diag(A)]/N;

end
