function [ Ec, Ev ] = BandD( N,Phi,Xi,Eg,Vr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [ c ] = constants();
    Ec = zeros(N+1,1);
    Ev = zeros(N+1,1);
    for i=1:N+1
        Ec(i) = -Xi(i)-Phi(i)*c.Vt+Vr;
        Ev(i) = Ec(i) - Eg(i);
    end   
    Ec = Ec - Ev(1);
    Ev = Ev - Ev(1);

end

