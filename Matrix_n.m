function [ SE ] = Matrix_n( N,Phi,p,h,Theta_n,Pbi_n,Pbi_p,mun,B,SRH )
% Build matrix for n
%   Detailed explanation goes here
% Massives defenition
    E1 = zeros((N-1)*3-2,1);
    E2 = zeros((N-1)*3-2,1);
    E = zeros((N-1)*3-2,1);
      
    i = 1;   
    E1(i) = i;
    E1(i+1) = i;
    E2(i) = i;
    E2(i+1) = i+1;
    E(i) =  -2/(h(i+1)+h(i))* (1/h(i)* mun(i+1)*BER(-Pbi_n+Phi(i))...
            +1/h(i+1)* mun(i+2)*BER(-Phi(i+1)+Phi(i))) + (B-SRH(i))*p(i);
    E(i+1) = 2/(h(i+1)+h(i))*1/h(i+1)* mun(i+2)*BER(-Phi(i)+Phi(i+1));
    for i=2:1:N-2
        E1((3*i-2)-1) = i;
        E1((3*i-2)) = i;
        E1((3*i-2)+1) = i;

        E2((3*i-2)-1) = i-1;
        E2((3*i-2)) = i;
        E2((3*i-2)+1) = i+1;

        E((3*i-2)-1) = 2/(h(i+1)+h(i))*1/h(i)* mun(i+1)*BER(-Phi(i)-Theta_n(i)+Phi(i-1)+Theta_n(i-1));
        E((3*i-2)) = -2/(h(i+1)+h(i))* (1/h(i)* mun(i+1)*BER(-Phi(i-1)-Theta_n(i-1)+Phi(i)+Theta_n(i))...
                     +1/h(i+1)* mun(i+2)*BER(-Phi(i+1)-Theta_n(i+1)+Phi(i)+Theta_n(i))) + (B-SRH(i))*p(i);
        E((3*i-2)+1) = 2/(h(i+1)+h(i))*1/h(i+1)* mun(i+2)*BER(-Phi(i)-Theta_n(i)+Phi(i+1)+Theta_n(i+1));
    end
    i = N-1;
    E1((3*i-2)) = i;
    E1((3*i-2)-1) = i;
    E2((3*i-2)) = i;
    E2((3*i-2)-1) = i-1;
    E((3*i-2)-1) = 2/(h(i+1)+h(i))*1/h(i)* mun(i+1)*BER(-Phi(i)+Phi(i-1));
    E((3*i-2)) = -2/(h(i+1)+h(i))* (1/h(i)* mun(i+1)*BER(-Phi(i-1)+Phi(i))...
                 +1/h(i+1)* mun(i+2)*BER(-Pbi_p+Phi(i))) + (B-SRH(i))*p(i);
    
    SE = sparse(E1,E2,E,N-1,N-1);

end

