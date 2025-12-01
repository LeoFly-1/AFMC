function [mseq]=mSequence(coef,registers)
%--------------------------------------------------------------------------
% By Yufei Li ---- 24.11.2025
% This function is used to generate m-sequecne
% coef: Polynomial coefficients, sorted in ascending order, column vector
% registers: Initial state of register, column vector
% Copyright Notice: If you refer to this code, please cite this paper
%***************************************************
coef=coef(2:end);
m=length(coef);
len=2^m-1;
backQ=0;
mseq=zeros(1,len);
for i=1:len
    mseq(i)=registers(m);
    backQ = mod(sum(coef.*registers) , 2);
    registers(2:length(registers)) = registers(1:length(registers)-1);
    registers(1)=backQ;
end
mseq=mseq';
end