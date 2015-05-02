function [ c ] = constants()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    c.kb = 1.38e-23;
    c.hb = 1.053e-34;
    c.m0 = 9.109e-31;
    c.T = 300;
    c.q = 1.6e-19;
    c.eps = 12.847;
    c.eps0 = 8.85e-14;
    c.Vt = c.kb*c.T/c.q;                  % potential normilize
    c.ni = 5e17;
    c.Ld = sqrt(c.eps0*c.Vt/c.q/c.ni);  % Length normilize
    c.gamma = 0.27;
    c.x = 0.2;
    c.AlGaAs = 0.8;
    c.P0 = -c.q*c.ni/c.Ld;
end

