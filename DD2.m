 clc;close all;
N = 350;    
UN = [0.537,0.54,0.545,0.55,0.5525,0.555,0.5575,0.56,0.5625,0.5650,0.5675,0.57,0.58,0.59];
Nn2 = length(UN)+1;
[ qd ] = QD_parameters();
[ c ] = constants();
[ l,lr ] = structure_parameters();
[ Nc,Nv,ni1,eps,P,Eg,taun,taup,mun,mup,Xi,C(:,1),Theta_n1,Theta_p1,X,h,b,L,D1 ] = StructureProfile( N );
bon = Boundary(D1,L,Theta_p1(N+1)*c.Vt);

%% TMP
ni = ones(N+1,1)*L(1).ni;
P = zeros(N+1,1);
% B = -2e-11*c.ni/c.Vt*c.Ld^2;
B = 0;
D_POISSON = 5e-4;%4
D_CONT = 5e-4;%5
%%
    Nqd = fix(qd.l/h(1));   % width of QD layer  fix(qd.l/h(1))
    Nqd2 = N/2-fix(Nqd/2);  % place for QD layer, in a middle
    Vce = qd.gen/Nc(Nqd2)/fermi(1/2,-(qd.dEc1-qd.en2)/c.Vt);
    Vcp = qd.gep/Nv(Nqd2)/fermi(1/2,-(qd.dEv1-qd.ep)/c.Vt);
%% Profile   
h = h/c.Ld;
t = 1;
Tol_n = 1e-7;
Tol_p = 1e-7;
mod_n = 2;
mod_p = 2;    
C1 = C;
jj = 1;
while(t < Nn2)
    Pbi_n = bon.Pbi_n - UN(t)/c.Vt; 
    Pbi_p = bon.Pbi_p + UN(t)/c.Vt;   
    Phi_b_n = -UN(t)/c.Vt;
    Phi_b_p = +UN(t)/c.Vt;
    [ Phi,Phi_n,Phi_p,n,p] = Initial(N,Pbi_n,Pbi_p,Phi_b_n,Phi_b_p,X,h,b,bon);
    Phi_n(N+1) = Phi_b_p;
    Phi_p(1) = Phi_b_n;

    G_iter = 0;    
    q = -0.5;
    J1 = 100;
    tmp_J = 1;
    TT = D_CONT;
    while((tmp_J > 1e-2) && (G_iter <1000))       %tmp > 1e7  
        % Iterative procedure POISSSON SOLVER  
        Phi_tmp = Phi;
        n_tmp = n;
        p_tmp = p;     
        tmp = 1;    
        iter = 0;
        Phi11 = Phi;   
        % all zeros
         while ((tmp > 1e-10) && (iter < 1000))
            [rho,d_rho,gPhi,n,p] = rho_f( N,h*c.Ld,c.eps0*eps(2:N+1),Phi*c.Vt,Phi_n,Phi_p,ni,C,Theta_n1,Theta_p1,bon.Vr );
            [ S,f ] = Matrix_Phi_f( N,h*c.Ld,c.q*rho(2:N),-c.q*d_rho(2:N),gPhi,P,c.eps0*eps(2:N+1) );
            Phi1 = S\(-f);
            tmp = max(abs(Phi1));
            Phi(2:N) = Phi(2:N) + D_POISSON*Phi1/c.Vt;      
            iter = iter+1;
            [q,R,fn2,fn1,fp,P1,P2] = F_QD_A2(Vce*n(N/2),Vcp*p(N/2));            
            C = C1;
            for i=1:Nqd
                C(Nqd2+i) = C(Nqd2+i) + qd.Ns/qd.l*(q-qd.P_A);  % because q<0 for -
            end               
        end                  
        if(R < 0)
            R = 0;
        end
        n = n/c.ni;
        p = p/c.ni;  
        [ n,p ] = CONTINUETY_SOLVER( N,X,Phi,n,p,Theta_n1,Theta_p1,mun,mup,ni1/c.ni,taun,taup,B,h,D_CONT,R,Nqd,Nqd2,Tol_n,Tol_p,mod_n,mod_p );
        % Let's calc Phi_n,Phi_p
        Phi_n = (Phi*c.Vt + Theta_n1*c.Vt - c.Vt*log(n./ni*c.ni)-bon.Vr)/c.Vt;
        Phi_p = (Phi*c.Vt - Theta_p1*c.Vt + c.Vt*log(p./ni*c.ni)-bon.Vr)/c.Vt;
        %% Break
        G_iter = G_iter + 1;
        [ Tol,J,Jn,Jp ] = CURRENT( N,Phi,n,p,h,Theta_n1,Theta_p1,mun,mup );
        tmp_J = abs(J-J1);
        J1 = J;
        figure(21)
        plot(X(1:N),c.Vt*c.ni/c.Ld*mun(1:N).*Jn,X(1:N),-c.Vt*c.ni/c.Ld*mup(1:N).*Jp)
    end
    hFactor = Vcp*p(N/2)/Vce/n(N/2);
    PhiPhi(:,t) = Phi;
    nn(:,t) = n;
    pp(:,t) = p;
    I(t) = J;
    [q,J,P1,P2,hFactor,Vce*n(N/2)*c.ni,Vcp*p(N/2)*c.ni]
    t = t+1;
end
