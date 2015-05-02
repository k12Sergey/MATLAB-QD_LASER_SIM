function [ n,p ] = CONTINUETY_SOLVER( N,X,Phi,n,p,Theta_n1,Theta_p1,mun,mup,ni,tau_n,tau_p,B,h,D,R,Nqd,Nqd2,Tol_n,Tol_p,mod_n,mod_p )
% We solve continuety equation here
%   Detailed explanation goes here
    [ c ] = constants();
    [ qd ] = QD_parameters();
    n11 = n;
    p11 = p;
    Pbi_n = Phi(1);
    Pbi_p = Phi(N+1);
    tmp = 1;
    iter = 0;
%     SRH = 1./(tau_n.*(p11+ni)+tau_p.*(n11+ni))/c.Vt*c.Ld^2;
    SRH = zeros(N+1,1);
%     figure(21)
    while ((tmp > Tol_n) && (iter < 100))% 1e-7 100
        [ F ] = C_Newton( N,h,Phi,n,p11,mup,mun,Theta_n1,Theta_p1,ni/c.ni,B,SRH,1 );
        %% QD TEMP
        for i=1:Nqd
            F(Nqd2+i) = F(Nqd2+i) - (B-SRH(i+1))*n(Nqd2+i)*p11(Nqd2+i)-SRH(Nqd2+i)*ni(Nqd2+i)^2 - R*c.Ld^2/c.Vt*qd.Ns/qd.l/c.ni;
        end           
        if(rem(iter,mod_n) == 0)
            [ SE ] = Matrix_n( N,Phi(2:N),p11(2:N),h,Theta_n1(2:N),Pbi_n,Pbi_p,mun,B,SRH(2:N) );
            for i=1:Nqd%                     
                SE(Nqd2+i,Nqd2+i) = SE(Nqd2+i,Nqd2+i) - (B-SRH(i+1))*p11(Nqd2+i); 
            end                
        end
        n1 = -SE\F;
        tmp = max(abs(n1));
        n(2:N) = n(2:N) + D*n1;
        iter = iter+1;
%         plot(X,log10(n))
    end 
    %% CONTINUETY HOLES
    tmp = 1;
    iter = 0;
    while ((tmp > Tol_p) && (iter < 100))% 1e-7 300
        [ F ] = C_Newton( N,h,Phi,n11,p,mup,mun,Theta_n1,Theta_p1,ni/c.ni,B,SRH,2 );
        %% QD TEMP
        for i=1:Nqd
            F(Nqd2+i) = F(Nqd2+i) - (B-SRH(i+1))*n11(Nqd2+i)*p(Nqd2+i)-SRH(Nqd2+i)*ni(Nqd2+i)^2 - R*c.Ld^2/c.Vt*qd.Ns/qd.l/c.ni; 
        end  
        if(rem(iter,mod_p) == 0)
            [ SH ] = Matrix_p( N,Phi(2:N),n11(2:N),h,Theta_p1(2:N),Pbi_n,Pbi_p,mup,B,SRH(2:N) );
            for i=1:Nqd
                SH(Nqd2+i,Nqd2+i) = SH(Nqd2+i,Nqd2+i) - (B-SRH(i+1))*n11(Nqd2+i);
            end
        end
        p1 = -SH\F;  
        tmp = max(abs(p1));
        p(2:N) = p(2:N) + D*p1;  
        iter = iter+1;  
    end 

end

