function [ Tol,J,Jn,Jp ] = CURRENT( N,Phi,n,p,h,Theta_n,Theta_p,mun,mup )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    [ c ] = constants();
    Jn = zeros(N,1);
    Jp = zeros(N,1);
    for i = 1:N
        Jn(i) = c.q/h(i)*(BER((-Phi(i+1)-Theta_n(i+1)+Phi(i)+Theta_n(i)))*n(i)-BER((Phi(i+1)+Theta_n(i+1)-Phi(i)-Theta_n(i)))*n(i+1));
        Jp(i) = c.q/h(i)*(BER((Phi(i+1)-Theta_p(i+1)-Phi(i)+Theta_p(i)))*p(i)-BER((-Phi(i+1)+Theta_p(i+1)+Phi(i)-Theta_p(i)))*p(i+1));
    end
    J_FULL = c.Vt*c.ni/c.Ld*mun(1:N).*Jn-c.Vt*c.ni/c.Ld*mup(1:N).*Jp;
    tmp_J = J_FULL - J_FULL(N/2);
    Tol = max(abs(tmp_J))/abs(J_FULL(N/2))*100;
    J = J_FULL(N/2);
end

