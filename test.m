clc;clear all;close all;
%% ========================================================================
% By Yufei Li, 24.11.2025
% This is a test demo for evaluating AFMC
% Copyright Notice: If you refer to this code, please cite this paper
%% ========================================================================
n = 3;                              % order of m-sequence
coef = [1 1 0 1].';                 % coefficient, ascending order
reg = [1 0 0].';                    % initial state
mSeq = GenerateMSequence(coef,reg); % generate the m-sequence

M_matrix = mSeq;
for i = 2:length(mSeq)
    M_matrix = [M_matrix circshift(M_matrix(:,i-1),-1)];
end
B = ones(2^n,2^n);
B(2:end,2:end) = 1-2*M_matrix;      % generate the bipolar augmented m-sequence matrix B

[Pr,Pc] = GenerateP(mSeq,coef);
H = hadamard(2^n);
B1 = Pr*H*Pc;                       % evaluate the permutation

isequal(B,B1)
%% ========================================================================