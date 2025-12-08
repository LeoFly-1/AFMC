function [Pr,Pc] = GenerateP(mSeq,coef)
% -------------------------------------------------------------------------
% Copyright Notice: If you refer to this code, please cite this paper
% By Yufei Li, 24.11.2025
% This function is used to generate the permutation matrices.
% Input: mSeq, column vector
% Output: row permutation matrix Pr and column permutation matrix Pc.
n = length(coef)-1;
S = mSeq;
for i = 2:n
    S(:,i) = circshift(S(:,i-1),-1);
end
cols = bi2de(S,'left-msb');
Pc = eye(2^n);
Pc = Pc(:,[1;cols+1]);

coef = flipud(coef(2:end));             % 
A = [zeros(n-1,1), eye(n-1); coef.'];   % A is the companion matrix of the coef
u = [1,zeros(1,n-1)];                   % u row vector, 1 x n
U = u;
iters = 2^n-1;
temp = u;
for i = 2:iters
    temp = mod(temp*A,2);
    U(i,:)=temp;
end
rows = bi2de(U,'left-msb');
Pr = eye(2^n);
Pr = Pr([1;rows+1],:);
end