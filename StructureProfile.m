function [ Nc,Nv,ni,eps1,P,Eg,taun,taup,mun,mup,Xi,C,Theta_n,Theta_p,X,h,b,L,D ] = StructureProfile( N )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
[ l,lr ] = structure_parameters();
[ c ] = constants();
% Pe = Piezo(N,X);
ni = zeros(N+1,1);
Nc = zeros(N+1,1);
Nv = zeros(N+1,1);
eps = zeros(N+1,1);
P = zeros(N+1,1);
Pe = zeros(N+1,1);
Eg = zeros(N+1,1);
taun = zeros(N+1,1);
taup = zeros(N+1,1);
mun = zeros(N+1,1);
mup = zeros(N+1,1);
Xi = zeros(N+1,1);
C = zeros(N+1,1);
Theta_n = zeros(N+1,1);
Theta_p = zeros(N+1,1);
%%STRUCTURE
%  b = [100,3,10,3,10,3,10,3,10,3,10,10,100]*1e-7; % nm
%  L = [l(1),l(2),l(1),l(2),l(1),l(2),l(1),l(2),l(1),l(2),l(1),l(3),l(1)];
%  D = [2,0,0.5,0,0.5,0,0.5,0,0.5,0,-2,-50,-2]*1e17;
b = [200,300,200.01]*1e-7; % nm
L = [l(6),l(7),l(6)];
D = [0.5,0,-0.5]*1e18;
 
[ X,h ] = grid( sum(b),N );
for i=1:N+1
    %% Choose the crierion
    temp = b(1);
    j = 1;
    while(true)
        if(X(i)<=temp)
            tmp_l = j;
            break;
        end
        j = j + 1;
        temp = temp + b(j);
    end
    %% Fill
    Nc(i) = L(tmp_l).Nc;
    Nv(i) = L(tmp_l).Nv;
    ni(i) = L(tmp_l).ni;
    eps(i) = L(tmp_l).eps;
    Pe(i) = 2*(L(tmp_l).a-l(lr).a)/L(tmp_l).a*(L(tmp_l).e31-L(tmp_l).e33*L(tmp_l).C13/L(tmp_l).C33);
    P(i) = L(tmp_l).P0 - Pe(i);
    Eg(i) = L(tmp_l).Eg;
    taun(i) = L(tmp_l).tau_n;
    taup(i) = L(tmp_l).tau_p;
    mun(i) = L(tmp_l).mun;
    mup(i) = L(tmp_l).mup;
    Xi(i) = L(tmp_l).Xi;
    C(i) = D(tmp_l); 
    Theta_n(i) = (L(tmp_l).Xi-l(lr).Xi)+c.Vt*log(L(tmp_l).Nc/l(lr).Nc);% - it was +
    Theta_p(i) = (-(L(tmp_l).Xi-l(lr).Xi)-(L(tmp_l).Eg-l(lr).Eg))+c.Vt*log(L(tmp_l).Nv/l(lr).Nv);
end
Theta_n = Theta_n/c.Vt;
Theta_p = Theta_p/c.Vt;
% Finite Volumes    
for i=2:N+1
    eps1(i) = (eps(i-1)+eps(i))/2;            
end    
end


