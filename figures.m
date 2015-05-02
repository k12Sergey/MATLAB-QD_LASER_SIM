function [  ] = figures( N,Nn2,C,n,p,Phi,U,X,Jn,Jp )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    [ c ] = constants();
    [ l1,l2,l3,het,b ] = structure_parameters();

    [ Ec,Ev ] = BandD( N,Phi,X );
%% 
    figure(100)
    plot(X,Phi*c.Vt-het.Vr);
    xlabel('X, cm')
    ylabel('\Phi, V')
%%
    figure(200)
    plot(X,log10(n*c.ni),X,log10(p*c.ni));
    xlabel('X, cm')
    ylabel('log10(n),log10(p), cm-3')
    legend('n','p');
%%
    figure(300)
    plot(X,Ec-Ev(1),X,Ev-Ev(1));
    xlabel('X, cm')
    ylabel('E, eV')
    legend('Ec','Ev');
%%
    figure(400)
    plot(X,C,X,n,X,-p);
    xlabel('X, cm')
    ylabel('C,n,p, cm-3')
    legend('C','n','p');
%%
    figure(500)
    plot(X(1:N),Jn,X(1:N),Jp);
    xlabel('X, cm')
    ylabel('Jn,Jp  A/cm^2')
    legend('Jn','Jp');
end

