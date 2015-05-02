function [ R, n ] = Potential(F,kb,T,q,eps,eps0,Nd,Nc)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    a = Nd;
    b = Nc*kb*T;
    F1 = F*sqrt(eps*eps0/2);
%     E0 = fminbnd(@tmp_zero,0,-q,[],a,b,kb*T);
%     x1 = linspace(0,q,100);
    E0 = -kb*T*log(Nd/Nc);
    C = (Nd*E0 + Nc*kb*T*exp(-E0/kb/T));
    R = fsolve(@tmp_Potential,E0,[],C,F1,a,b,kb*T);
    n = Nc*exp(-R/T/kb);
%     R = R/q;    
end
function y = tmp_Potential(x,C,F,a,b,T)
    y = F^2-a*x-b.*exp(-x/T)+C;
end
% function y = tmp_min(x,a,b,T)
%     y = -a*x-b*exp(-x/T);
% end