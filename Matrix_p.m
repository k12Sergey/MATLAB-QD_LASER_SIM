function [ SH ] = Matrix_p( N,Phi,n,h,Theta_p,Pbi_n,Pbi_p,mup,B,SRH )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
% LEt it be like this mun(i-1/2) = mun(i), mun(i+1/2) = mun(i+1)
       
    H1 = zeros((N-1)*3-2,1);
    H2 = zeros((N-1)*3-2,1);
    H = zeros((N-1)*3-2,1);
    
    i = 1;   
    H1(i) = i;
    H1(i+1) = i;
    H2(i) = i;
    H2(i+1) = i+1;
    H(i) = -2/(h(i+1)+h(i))* (1/h(i)* mup(i+1)*BER(Pbi_n-Phi(i))...
           +1/h(i+1) *mup(i+2)*BER(Phi(i+1)-Phi(i))) + (B-SRH(i))*n(i);    
    H(i+1) = 2/(h(i+1)+h(i))*1/h(i+1)* mup(i+2)*BER(Phi(i)-Phi(i+1));
    for i=2:1:N-2
        H1((3*i-2)-1) = i;
        H1((3*i-2)) = i;
        H1((3*i-2)+1) = i;

        H2((3*i-2)-1) = i-1;
        H2((3*i-2)) = i;
        H2((3*i-2)+1) = i+1;

        H((3*i-2)-1) = 2/(h(i+1)+h(i))*1/h(i)* mup(i+1)*BER(Phi(i)-Theta_p(i)-Phi(i-1)+Theta_p(i-1));
        H((3*i-2)) = -2/(h(i+1)+h(i))* (1/h(i)* mup(i+1)*BER(Phi(i-1)-Theta_p(i-1)-Phi(i)+Theta_p(i)) ...
                     +1/h(i+1)* mup(i+2)*BER(Phi(i+1)-Theta_p(i+1)-Phi(i)+Theta_p(i))) + (B-SRH(i))*n(i);
        H((3*i-2)+1) = 2/(h(i+1)+h(i))*1/h(i+1)* mup(i+2)*BER(Phi(i)-Theta_p(i)-Phi(i+1)+Theta_p(i+1));
    end
    i = N-1;
    H1((3*i-2)) = i;
    H1((3*i-2)-1) = i;
    H2((3*i-2)) = i;
    H2((3*i-2)-1) = i-1;
    H((3*i-2)-1) = 2/(h(i+1)+h(i))*1/h(i)* mup(i+1)*BER(Phi(i)-Phi(i-1));
    H((3*i-2)) = -2/(h(i+1)+h(i))* (1/h(i)* mup(i+1)*BER(Phi(i-1)-Phi(i))...
                 +1/h(i+1)* mup(i+2)*BER(Pbi_p-Phi(i))) + (B-SRH(i))*n(i);             % I change the sign

    SH = sparse(H1,H2,H,N-1,N-1);

end

