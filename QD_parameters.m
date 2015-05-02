function [ qd ] = QD_parameters()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    qd.gen = 219e9;
    qd.gep = 19e9;
    qd.gr1 = 1e9;
    qd.gr2 = 1e9;
    qd.ge12 = 125e9;
    qd.ge21 = 1250e9;
    qd.gamma1 = 0.6;% 0.55
    qd.gamma2 = 0.333;% 0.22 
    qd.D1 = 2;
    qd.D2 = 4;
    qd.dEc1 = 0.31;
    qd.dEv1 = 0.31;
    qd.en1 = 0.09;
    qd.en2 = 0.15;
    qd.ep = 0.09;
    qd.l = 10e-7; % 10 ncm
    %%
    qd.Ns = 10e10; % 1/cm^2
    qd.P_A = 0;
end

