function [ E,gE ] = Tests( N,X,h,Phi,n,p,C,eps,mun,mup,Theta_n,Theta_p,B)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
	[ c ] = constants();
    %% Poisson
%     for i=2:N
%         E(i,:) = -(Phi(i)-Phi(i-1))/h(i)*c.Vt/c.Ld;
%     end
%     for i=3:N-1
%         gE(i,:) = (eps(i)*E(i)-eps(i-1)*E(i-1))/h(i)/c.Ld;        
%     end
%     figure(1000)
%     plot(X(3:N-1),gE(3:N-1),X(3:N-1),c.q/c.eps0*(p(2:N-2)*c.ni-n(2:N-2)*c.ni+C(2:N-2)));     
%     xlabel('X,cm')
%     ylabel('a.u')
%     legend('d/dx(eps d/dx\phi)','q*(p-n+C)')
%     title ('Poisson equation - test')
    %% Current
    for i = 1:N
        Jn(i) = c.q/h(1)*(BER((-Phi(i+1)-Theta_n(i+1)+Phi(i)+Theta_n(i)))*n(i)-BER((Phi(i+1)+Theta_n(i+1)-Phi(i)-Theta_n(i)))*n(i+1));
        Jp(i) = c.q/h(1)*(BER((Phi(i+1)-Theta_p(i+1)-Phi(i)+Theta_p(i)))*p(i)-BER((-Phi(i+1)+Theta_p(i+1)+Phi(i)-Theta_p(i)))*p(i+1));
    end
%     for i = 2:N
%         dJn(i,:) = (Jn(i)-Jn(i-1))/h(i);
%         dJp(i,:) = (Jp(i)-Jp(i-1))/h(i);
%         F(i) = B*n(i)*p(i); %*c.ni^2*c.q*c.Vt/c.ni/c.Ld
%     end
    figure(2000)
    plot(X(1:N),c.Vt*c.ni./c.Ld*mun(2:N).*Jn,X(1:N),-c.Vt*c.ni./c.Ld*mup(2:N).*Jp);
%     xlabel('X,cm')
%     ylabel('J, A/cm^2')
%     legend('J_n','J_p')
%     figure(3000)
%     plot(X(2:N),mun(10)*dJn(2:N)/c.q,X(2:N),F(2:N)) 
%     xlabel('X,cm')
%     legend('div J_n','Bnp')
%     title('current equation, arbitrary units - test')
    figure(4000)
    plot(X(1:N),c.Vt*c.ni/c.Ld*mun(2:N)*Jn-c.Vt*c.ni/c.Ld*mup(2:N)*Jp)
end

