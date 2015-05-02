function [ q,R,fe2,fe1,fh,P1,P2 ] = F_QD_A2(gcn,gcp)
[ qd ] = QD_parameters();
Pe = gcn/(gcn+qd.gen);
Ph = gcp/(gcp+qd.gep);
D = qd.D2/qd.D1;
C = 1/(1+D^(-1))*(gcn+qd.gen)/(gcp+qd.gep);
eta = (gcn+qd.gen)/qd.ge21;
xi = qd.ge12/qd.ge21;
g1 = qd.gamma1;
g2 = qd.gamma2;
%% QD::NL    
    P1 = 0;
    P2 = 0;
    fe2 = Pe;
    fe1 = Pe/((1-xi)*Pe+xi);
    fh = Ph;        
    if(fe1+fh-1 > qd.gamma1 || fe2+fh-1 > qd.gamma2)
        P1=-(1/(2*C*(-1 + xi)))*(eta-g1+C*Pe+Ph+xi+C*xi+g1*xi-C*Pe*xi - ...
            Ph*xi-sqrt(4*(1+g1+C*Pe+eta*(1+g1-Ph)-Ph)*(-1+xi)+(2+eta+g1+C*Pe-Ph-xi+C*xi-g1*xi- ... 
            C*Pe*xi+Ph*xi)^2));
        P2 = 0;
        fe2 = Pe-P1;
        fe1 = ((1+eta)*P1-Pe)/((xi-1)*(Pe-P1)-xi);
        fh = Ph-C*P1;  
        if(fe2+fh-1 > qd.gamma2)
            fe2 = (C*Pe-Ph+1+g2)/(C+1);
            fe1 = fe2/(fe2*(1-xi)+xi);
            fh = 1+qd.gamma2-fe2;
            P2 = (Ph-fh)/C;
            if(fe1+fh-1 > qd.gamma1)
                fh = (Ph-C*(Pe-1-qd.gamma2))/(C+1);
                fe2 = 1+qd.gamma2 - fh;
                fe1 = 1+qd.gamma1 - fh;
                P1 = (fe2*(1-fe1)-xi*fe1*(1-fe2))/eta;
                P2 = (Ph-fh)/C-P1;
            else
                P1 = 0;
            end
        end
    end                
% %     % test
%     P1 = P1*(gcn+qd.gen)*qd.D2;
%     P2 = P2*(gcn+qd.gen)*qd.D2;
%     dd = gcn*(1-fe2)-qd.gen*fe2-qd.ge21*fe2*(1-fe1)+qd.ge12*fe1*(1-fe2)-P2/qd.D2
%     dd1 = qd.D2/qd.D1*(qd.ge21*fe2*(1-fe1)-qd.ge12*fe1*(1-fe2))-P1/qd.D1
%     dd2 = gcp*(1-fh)-qd.gep*fh-(P1+P2)/6
    %output
    P1 = P1*(gcn+qd.gen)*qd.D2*4e-3*qd.Ns*983*1e-3*1.6e-19;
    P2 = P2*(gcn+qd.gen)*qd.D2*4e-3*qd.Ns*1050*1e-3*1.6e-19;
    q = (qd.D1+qd.D2)*fh-qd.D1*fe1-qd.D2*fe2;    
    R = 6*(gcp*(1-fh)-qd.gep*fh)+2*qd.gr1*fe1*fe2+4*qd.gr2*fe2*fe1;
    
end
