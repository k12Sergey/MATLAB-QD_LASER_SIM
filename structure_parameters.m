function [ l, l_r] = structure_parameters()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Constants
    [ c ] = constants();
% band structure parameters
    %% layer 1
    l1.Eg = 3.51;
    l1.Xi = 1.96;
    l1.mn = 0.2*c.m0;
    l1.mlh = 0.3*c.m0;
    l1.mhh = 0.02*c.m0;
    l1.Nc = 2*(l1.mn*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l1.Nv = 2*((l1.mlh^(3/2)+l1.mhh^(3/2))^(2/3)*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l1.ni = sqrt(l1.Nc*l1.Nv)*exp(-l1.Eg/2/c.Vt);
    l1.mun = 7000;
    l1.mup = 400;
    l1.tau_n = 10e-9;
    l1.tau_p = 100e-9; 
    l1.eps = 8.9;   
    l1.Ed = 0.013; 
    l1.P0 = -0.029*1e-4;
    l1.a =  0.3188*1e-7; % along z axis ?!
    l1.e31 = -0.33*1e-4;
    l1.e33 = 0.65*1e-4;
    l1.C13 = 105;
    l1.C33 = 395;
    l1.B = 2.4e-11;
    l1.Cn = 0;
    l1.Cp = 0;
   %% l2 In0.2Ga0.8N
    l2.Eg = c.x*0.69+(1-c.x)*3.51-1.2*c.x*(1-c.x);
    l2.Xi = c.x*3.85+(1-c.x)*1.96+0.84*c.x*(1-c.x);
    l2.mn = 0.2*c.m0;
    l2.mlh = 0.01*c.m0;
    l2.mhh = 1.15*c.m0;
    l2.Nc = 2*(l2.mn*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l2.Nv = 2*((l2.mlh^(3/2)+l2.mhh^(3/2))^(2/3)*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l2.ni = sqrt(l2.Nc*l2.Nv)*exp(-l2.Eg/2/c.Vt);
    l2.mun = 7000;
    l2.mup = 400;
    l2.tau_n = 10e-9;
    l2.tau_p = 100e-6;    % 1e-6 
    l2.eps = 14.6;  %8.9
    l2.Ed = 0;
    l2.P0 = (c.x*(-0.032)+(1-c.x)*(-0.029))*1e-4;
    l2.a =  (c.x*0.354+(1-c.x)*0.3188)*1e-7; % ncm
    l2.e31 = (c.x*(-0.22)+(1-c.x)*(-0.33))*1e-4;
    l2.e33 = (c.x*(0.43)+(1-c.x)*(0.65))*1e-4;
    l2.C13 = c.x*(95)+(1-c.x)*(105);
    l2.C33 = c.x*(200)+(1-c.x)*(395);
    l2.B = c.x*2.4e-11+(1-c.x)*6.6e-12;
    l2.Cn = (1-c.x)*2.5e-30;
    l2.Cp = (1-c.x)*2.5e-30;
   %% l3 AlGaN
    c.x = 0.15;
    l3.Eg = c.x*6.25+(1-c.x)*3.51-1*c.x*(1-c.x);
    l3.Xi = c.x*0+(1-c.x)*1.96+0.7*c.x*(1-c.x);
    l3.mn = 0.2*c.m0;
    l3.mlh = 1.15*c.m0;
    l3.mhh = 0.2*c.m0;
    l3.Nc = 2*(l3.mn*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l3.Nv = 2*((l3.mlh^(3/2)+l3.mhh^(3/2))^(2/3)*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l3.ni = sqrt(l3.Nc*l3.Nv)*exp(-l3.Eg/2/c.Vt);
    l3.mun = 7000;
    l3.mup = 400;
    l3.tau_n = 10e-9;
    l3.tau_p = 100e-6;    % 1e-6 
    l3.eps = 8.5;  %8.9
    l3.Ed = 0;
    l3.P0 = (c.x*(-0.081)+(1-c.x)*(-0.029))*1e-4;
    l3.a =  (c.x*0.3112+(1-c.x)*0.3188)*1e-7; % ncm
    l3.e31 = (c.x*(-0.58)+(1-c.x)*(-0.33))*1e-4;
    l3.e33 = (c.x*(1.55)+(1-c.x)*(0.65))*1e-4;
    l3.C13 = c.x*(115)+(1-c.x)*(105);
    l3.C33 = c.x*(385)+(1-c.x)*(395);
    l3.B = c.x*2.4e-11+(1-c.x)*6.6e-12;
    l3.Cn = (1-c.x)*2.5e-30;
    l3.Cp = (1-c.x)*2.5e-30;
   %%  GaSb
    l4.Eg = 0.6;
    l4.Xi = 4.1;
    l4.mn = 0.041*c.m0;
    l4.mlh = 0.05*c.m0;
    l4.mhh = 0.4*c.m0;
    l4.Nc = 2*(l4.mn*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l4.Nv = 2*((l4.mlh^(3/2)+l4.mhh^(3/2))^(2/3)*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l4.ni = sqrt(l4.Nc*l4.Nv)*exp(-l4.Eg/2/c.Vt);
    l4.mun = 5000;
    l4.mup = 2000;
    l4.tau_n = 10e-9;
    l4.tau_p = 100e-9; 
    l4.eps = 15;   
    l4.Ed = 0.013; 
    l4.P0 = 0;
    l4.a =  0; % along z axis ?!
    l4.e31 = 0;
    l4.e33 = 0;
    l4.C13 = 0;
    l4.C33 = 0;
    l4.B = 0;
    l4.Cn = 0;
    l4.Cp = 0;
   %% InAs
    l5.Eg = 0.35;
    l5.Xi = 4.9;
    l5.mn = 0.023*c.m0;
    l5.mlh = 0.026*c.m0;
    l5.mhh = 0.41*c.m0;
    l5.Nc = 2*(l5.mn*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l5.Nv = 2*((l5.mlh^(3/2)+l5.mhh^(3/2))^(2/3)*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l5.ni = sqrt(l5.Nc*l5.Nv)*exp(-l5.Eg/2/c.Vt);
    l5.mun = 7000;
    l5.mup = 400;
    l5.tau_n = 10e-9;
    l5.tau_p = 100e-9;
    l5.eps = 15;
    l5.Ed = 0.17;
    l5.P0 = -0;
    l5.a =  0; %ncm
    l5.e31 = -0;
    l5.e33 = 0;
    l5.C13 = 0;
    l5.C33 = 0;
    l5.B = 0;
    l5.Cn = 0;
    l5.Cp = 0;
    %% AlGaAs
    l6.Eg = 1.9+0.125*c.AlGaAs+0.143*c.AlGaAs^2;
    l6.Xi = 3.64-0.14*c.AlGaAs;
    l6.mn = (0.85-0.14*c.AlGaAs)*c.m0;
    l6.mlh = 0.021*c.m0;
    l6.mhh = 1.1*c.m0;
    l6.Nc = 2*(l6.mn*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l6.Nv = 2*((l6.mlh^(3/2)+l6.mhh^(3/2))^(2/3)*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l6.ni = sqrt(l6.Nc*l6.Nv)*exp(-l6.Eg/2/c.Vt);
    l6.mun = 7000;
    l6.mup = 400;
    l6.tau_n = 10e-8;
    l6.tau_p = 10e-6;    % 1e-6 
    l6.eps = 12.90-2.84*c.AlGaAs;  %8.9
    l6.Ed = 0;
    l6.P0 = 0;
    l6.a =  0; % ncm
    l6.e31 = 0;
    l6.e33 = 0;
    l6.C13 = 0;
    l6.C33 = 0;
    l6.B = 0;
    l6.Cn = 0;
    l6.Cp = 0;
    %%
    l7.Eg = 1.424;
    l7.Xi = 4.07;
    l7.mn = 0.063*c.m0;
    l7.mlh = 0.04*c.m0;
    l7.mhh = 0.51*c.m0;
    l7.Nc = 2*(l7.mn*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l7.Nv = 2*((l7.mlh^(3/2)+l7.mhh^(3/2))^(2/3)*c.kb*c.T/2/pi/c.hb^2)^(3/2)*1e-6;
    l7.ni = sqrt(l7.Nc*l7.Nv)*exp(-l7.Eg/2/c.Vt);
    l7.mun = 7000;   %212
    l7.mup = 400;    %70
    l7.tau_n = 10e-9;
    l7.tau_p = 10e-6;    % 1e-6 
    l7.eps = 12.9;  %8.9
    l7.Ed = 0;
    l7.P0 = 0;
    l7.a =  (c.x*0.354+(1-c.x)*0.3188)*1e-7; % ncm
    l7.e31 = 0;
    l7.e33 = 0;
    l7.C13 = 0;
    l7.C33 = 0;
    l7.B = 0;
    l7.Cn = 0;
    l7.Cp = 0;
    %%
    l = [l1,l2,l3,l4,l5,l6,l7];
%% Heterojunction values
    l_r = 6; 
    
  
    
end

