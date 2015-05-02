function [ bon ] = Boundary(D,L,T1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [l,lr] = structure_parameters();
    [ c ] = constants();
    bon.Vr = 0;% c.Vt*log(l(lr).Nc/D(1))-2*l(lr).Xi-l(lr).Eg/2-c.Vt/2*log(l(lr).Nv/l(lr).Nc);
%     options = optimset('MaxFunEvals',400000,'MaxIter',100000,'TolFun',1e-16);
%     T01 = [0.049,l1.N/c.ni,l1.ni^2/l1.N/c.ni];
%     T1 = fsolve(@BC,T01,options,l1.Nc/c.ni,l1.Nv/c.ni,l1.N/c.ni,l1.Eg,l1.Ed);
%     T02 = [0,l3.ni^2/l3.N/c.ni,l3.N/c.ni];
%     T2 = fsolve(@BV,T02,options,l3.Nc/c.ni,l3.Nv/c.ni,l3.N/c.ni,l3.Eg,l3.Ea);
%     %% Boundary
%     bon.nD = T1(2);
%     bon.pD = T1(3);
%     bon.nA = T2(2);
%     bon.pA = T2(3);
%     bon.Pbi_n = (T1(1)-l1.Xi-het.Vr)/c.Vt;
%     bon.Pbi_p = (-T2(1)-l1.Xi-het.Vr-l3.Eg)/c.Vt;

%% Bolzman Ohmic condition    
    N = length(L);
    if((D(1)<0) && (D(N)>0))
        bon.pD = abs(D(1))/c.ni;
        bon.nD = (L(1).ni^2/abs(D(1)))/c.ni;
        bon.pA = (L(N).ni^2/D(N))/c.ni;
        bon.nA = D(N)/c.ni;
        bon.Pbi_n = (-c.Vt*log(abs(D(1))/L(1).ni)+bon.Vr)/c.Vt;   %(c.Vt*log(l1.N/l1.ni)+het.Vr)/c.Vt;    
        bon.Pbi_p = (+c.Vt*log(abs(D(N))/L(N).ni)+T1+bon.Vr)/c.Vt;    %(-c.Vt*log(l3.N/l3.ni)+het.Vr)/c.Vt;
    elseif((D(1)<0) && (D(N)<0))
        bon.pD = abs(D(1))/c.ni;
        bon.nD = (L(1).ni^2/abs(D(1)))/c.ni;
        bon.nA = (L(N).ni^2/abs(D(N)))/c.ni;
        bon.pA = abs(D(N))/c.ni;
        bon.Pbi_n = (-c.Vt*log(abs(D(1))/L(1).ni)+bon.Vr)/c.Vt;   %(c.Vt*log(l1.N/l1.ni)+het.Vr)/c.Vt;    
        bon.Pbi_p = (-c.Vt*log(abs(D(N))/L(1).ni)+T1+bon.Vr)/c.Vt;    %(-c.Vt*log(l3.N/l3.ni)+het.Vr)/c.Vt;
    elseif((D(1)>0) && (D(N)>0))
        bon.nD = abs(D(1))/c.ni;
        bon.pD = (L(1).ni^2/abs(D(1)))/c.ni;
        bon.pA = (L(N).ni^2/abs(D(N)))/c.ni;
        bon.nA = abs(D(N))/c.ni;
        bon.Pbi_n = (c.Vt*log(abs(D(1))/L(1).ni)+bon.Vr)/c.Vt;   %(c.Vt*log(l1.N/l1.ni)+het.Vr)/c.Vt;    
        bon.Pbi_p = (c.Vt*log(abs(D(N))/L(1).ni)-T1+bon.Vr)/c.Vt;    %(-c.Vt*log(l3.N/l3.ni)+het.Vr)/c.Vt;
    elseif((D(1)>0) && (D(N)<0))
        bon.nD = abs(D(1))/c.ni;
        bon.pD = (L(1).ni^2/abs(D(1)))/c.ni;
        bon.nA = (L(N).ni^2/abs(D(N)))/c.ni;
        bon.pA = abs(D(N))/c.ni;
        bon.Pbi_n = (c.Vt*log(abs(D(1))/L(1).ni)+bon.Vr)/c.Vt;   %(c.Vt*log(l1.N/l1.ni)+het.Vr)/c.Vt;    
        bon.Pbi_p = (-c.Vt*log(abs(D(N))/L(1).ni)+T1+bon.Vr)/c.Vt;    %(-c.Vt*log(l3.N/l3.ni)+het.Vr)/c.Vt;
    end
    bon1 = bon;   
end

function y = BC(x,Nc,Nv,Nd,Eg,Ed)
    [ c ] = constants();
    y = [-Nc*fermi(1/2,x(1)/c.Vt)+Nv*fermi(1/2,(-x(1)-Eg)/c.Vt)+Nd/(1+2*exp((-FF(x(2)/Nc,1)+Ed)/c.Vt)), x(2)-Nc*fermi(1/2,x(1)/c.Vt), x(3)-Nv*fermi(1/2,(-x(1)-Eg)/c.Vt)];%/(1+2*exp((-FF(x(2)/Nc,1)+Ed)/c.Vt))
end
function y = BV(x,Nc,Nv,Na,Eg,Ea)
    [ c ] = constants();
    y = [Nv*fermi(1/2,x(1)/c.Vt)-Nc*fermi(1/2,(-x(1)-Eg)/c.Vt)-Na/(1+4*exp((FF(x(3)/Nv,2)+Ea)/c.Vt)), x(2)-Nc*fermi(1/2,(-x(1)-Eg)/c.Vt), x(3)-Nv*fermi(1/2,x(1)/c.Vt)];%/(1+4*exp((FF(x(3)/Nv,2)+Ea)/c.Vt))
end


