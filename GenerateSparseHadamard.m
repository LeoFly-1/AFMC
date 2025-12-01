function HH = GenerateSparseHadamard(n)
% -------------------------------------------------------------------------
% By Yufei Li, 24.11.2025
% This function is used to generate Sparse matrices of Hadamard matrix.
% Input: n is the order of the Hadamard matrix H, H is 2^n x 2^n.
% Output: HH, is a cell array, has n sparse matrices.
idx = zeros(n,2^n);
N = 2^n;
for i=0:n-1
    A=reshape(0:N-1,2^(n-i),2^i);
    B=A;
    B(1:end/2,:)=A(1:2:end,:);
    B(end/2+1:end,:)=A(2:2:end,:);
    idx(i+1,:)=B(:)+1;
end

H0=(kron(speye(N/2),[1 1;1,-1]));
tempeye=speye(N);

for i=1:n
    temp=H0(idx(i,:),idx(i,:));
    if i==1
        temp=(temp*tempeye);
    end
    HH{i}=(temp);
end

end