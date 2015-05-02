function [ Phi,n,p,R,q,fn,fp ] = POISSON_SOLVER( N,Phi,Phi_n,Phi_p,Vce,Vcp,C,P,Theta_n1,Theta_p1,ni,eps,Vr,h,D,Nqd,Nqd2 )
% We solve Poisson equation in here
%   Detailed explanation goes here
    [ c ] = constants();
    [ qd ] = QD_parameters();
    C1 = C;
    tmp = 1;    
    iter = 0;
    Phi11 = Phi;
    while ((tmp > 1e-10) && (iter < 1000))
        [rho,d_rho,gPhi,n,p] = rho_f( N,h*c.Ld,c.eps0*eps(2:N+1),Phi*c.Vt,Phi_n,Phi_p,ni,C,Theta_n1,Theta_p1,Vr );
        [ S,f ] = Matrix_Phi_f( N,h*c.Ld,c.q*rho(2:N),-c.q*d_rho(2:N),gPhi,P,c.eps0*eps(2:N+1) );
        Phi1 = S\(-f);
        tmp = max(abs(D*Phi1));
        Phi(2:N) = Phi(2:N) + D*Phi1/c.Vt;      
        iter = iter+1;
        [q,R,fn,fp] = F_QD(Vce*n(N/2),Vcp*p(N/2));
        C = C1;
        for i=1:Nqd
            C(Nqd2+i) = C(Nqd2+i) + qd.Ns/qd.l*q;  % because q<0 for -
        end               
    end      
    n = n/c.ni;
    p = p/c.ni;
end

