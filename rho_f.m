function [ y,f,g,n,p ] = rho_f( N,h,eps,Phi,Phi_n,Phi_p,ni,C,Theta_n,Theta_p,Vr )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    y = Destribution( Phi,ni,2,1,Phi_n,Phi_p,Theta_n,Theta_p,Vr )-Destribution( Phi,ni,1,1,Phi_n,Phi_p,Theta_n,Theta_p,Vr )+C;
    f = Destribution( Phi,ni,2,2,Phi_n,Phi_p,Theta_n,Theta_p,Vr )-Destribution( Phi,ni,1,2,Phi_n,Phi_p,Theta_n,Theta_p,Vr );

    n = Destribution(Phi,ni,1,1,Phi_n,Phi_p,Theta_n,Theta_p,Vr );
    p = Destribution( Phi,ni,2,1,Phi_n,Phi_p,Theta_n,Theta_p,Vr );
    
    g = zeros(N-1,1);
    for i=1:N-1
        g(i) = (eps(i)/h(i)*Phi(i)-(eps(i+1)/h(i+1)...
                +eps(i)/h(i))*Phi(i+1)+eps(i+1)/h(i+1)*Phi(i+2));
    end
end

