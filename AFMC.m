function [Y] = AFMC(sig,F,Pr,Pc,Spc,n,path)
% -------------------------------------------------------------------------
% Copyright Notice: If you refer to this code, please cite this paper
% By Yufei Li, 24.11.2025
% This function is the AFMC algorithm.
% Input: sig, column vector, (2^n-1)*Spc.
%        F, frequency compresation mat, length(sig) x Ka
%        Pr,Pc, permutation matices, 2^n x 2^n
%        Spc, Samples per chip
%        n, order of m-sequence
%        path: 1 or 2, 1 is the fwht, 2 is the sparse matrix;
% Output: Y
N = 2^n;
Y = zeros(length(sig),size(F,2));
ck = zeros(N,1);
if isequal(path,1)
    for i = 1:size(F,2)
        rk = sig.*F(:,i);
        for j = 1:Spc
            ck(2:end) = (sum(reshape(rk,Spc,N-1))).';
            yk = Pr*fwht((Pc*ck),N,'hadamard');
            Y(j:Spc:end,i) = abs(yk);
            rk = circshift(rk,-1);
        end
    end
elseif isequal(path,2)
    HH = GenerateSparseHadamard(n);
    for i = 1:size(F,2)
        rk = sig.*F(:,i);
        for j = 1:Spc
            ck(2:end) = (sum(reshape(rk,Spc,N-1))).';
            yk = Pc*ck;
            for k = 1:n
                yk = HH{k}*yk;
            end
            yk = Pr*yk;
            Y(j:Spc:end,i) = abs(yk);
            rk = circshift(rk,-1);
        end
    end
end
end