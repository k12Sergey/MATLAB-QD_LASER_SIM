function [ S,f ] = Matrix_Phi_f( N,h,rho,d_rho,gPhi,PL,eps )
% Buid linear system and right side
%  Potential
    [ c ] = constants();
% Defenition
    P1 = zeros((N-1)*3-2,1);
    P2 = zeros((N-1)*3-2,1);
    P = zeros((N-1)*3-2,1);
    f = zeros(N-1,1);
% System
    i = 1;
    P1(i) = i;
    P1(i+1) = i;
    P2(i) = i;
    P2(i+1) = i+1;
    P(i) = -(eps(i+1)/h(i+1)+eps(i)/h(i)+d_rho(i)*h(i));
    P(i+1) = eps(i+1)/h(i+1);
    for i=2:1:N-2
        P1((3*i-2)-1) = i;
        P1((3*i-2)) = i;
        P1((3*i-2)+1) = i;

        P2((3*i-2)-1) = i-1;
        P2((3*i-2)) = i;
        P2((3*i-2)+1) = i+1;

        P((3*i-2)-1) = eps(i)/h(i);
        P((3*i-2)) = -(eps(i+1)/h(i+1)+eps(i)/h(i)+d_rho(i)*h(i));
        P((3*i-2)+1) = eps(i+1)/h(i+1);
    end
    i = N-1;
    P1((3*i-2)) = i;
    P1((3*i-2)-1) = i;
    P2((3*i-2)) = i;
    P2((3*i-2)-1) = i-1;
    P((3*i-2)-1) = eps(i)/h(i);
    P((3*i-2)) = -(eps(i+1)/h(i+1)+eps(i)/h(i)+d_rho(i)*h(i));
    
% Right side with boundaries i=1,N-1
    i = 1;
    f(i) = rho(i)*h(i)+gPhi(i)-(PL(i+1)-PL(i));%-Pbi_n*2/(h(i)+h(i+1))*eps(i)/h(i)
    for i = 2:1:N-2
        f(i) = rho(i)*h(i)+gPhi(i)-(PL(i+1)-PL(i));
    end
    i = N-1;
    f(N-1) = rho(i)*h(i)+gPhi(i)-(PL(i+1)-PL(i));%-Pbi_p*2/(h(i)+h(i+1))*eps(i+1)/h(i+1)   
    
    S = sparse(P1,P2,P,N-1,N-1);
    
end