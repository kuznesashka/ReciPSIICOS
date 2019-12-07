function iM = tihinv(M, alpha)

if(nargin==1)
    alpha = 0.00001;
end;

MeanDiag = mean(abs(diag(M)));

MTikh = M+MeanDiag*alpha*eye(size(M,1));

iM = inv(MTikh);