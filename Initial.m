function [ Phi,Phi_n,Phi_p,n,p ] = Initial(N,Pbi_n,Pbi_p,Phi_b_n,Phi_b_p,X,h,b,bon)
% Initial
%   Detailed explanation goes here
%     [ l,lr ] = structure_parameters();
    [ c ] = constants();
%     [ X,h ] = grid(sum(b),N );
%     X = X*1e-4;
%     h = h*1e-4;
%   
% Initial potential and doping profile
    Phi = zeros(N,1);   % potential
    Phi_n = zeros(N,1); % quasi-Fermi electrons
    Phi_p = zeros(N,1); % quasi-Fermi holes
    n = zeros(N,1);     % electron concentratio
    p = zeros(N,1);     % hole concentration

    for i=1:N+1
        if(X(i) <= b(1))
           Phi(i) = Pbi_n;  
           Phi_n(i) = Phi_b_n;
           Phi_p(i) = Phi_b_p;
        elseif(X(i) <= sum(b(1:length(b)-1)))
           Phi(i) = (Pbi_n + Pbi_p)/2;
           Phi_n(i) = Phi_b_n;
           Phi_p(i) = Phi_b_p;
        elseif(X(i) <= sum(b))
           Phi(i) = Pbi_p;
           Phi_n(i) = Phi_b_n;
           Phi_p(i) = Phi_b_p;
        end
    end
   
% INitial concentrations
    for i = 1:N+1
        if(X(i) <= b(1)) 
            n(i) = bon.nD;
            p(i) = bon.pD;
        elseif(X(i) <= sum(b(1:length(b)-1)))
            n(i) = 0/c.ni;
            p(i) = 0/c.ni;
        elseif(X(i) <= sum(b))
            n(i) = bon.nA;
            p(i) = bon.pA;
        end
    end
end

