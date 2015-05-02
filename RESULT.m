N = 350;
[ qd ] = QD_parameters();
[ c ] = constants();
[ Nc,Nv,ni1,eps,P,Eg,taun,taup,mun,mup,Xi,C(:,1),Theta_n1,Theta_p1,X,h,b,L,D1 ] = StructureProfile( N );
Nqd = fix(qd.l/h(1));   % width of QD layer  fix(qd.l/h(1))
Nqd2 = N/2-fix(Nqd/2);  % place for QD layer, in a middle
Vce = qd.gen/Nc(Nqd2)/fermi(1/2,-(qd.dEc1-qd.en2)/c.Vt);
Vcp = qd.gep/Nv(Nqd2)/fermi(1/2,-(qd.dEv1-qd.ep)/c.Vt);
load('11.mat')
nn1 = nn;
pp1 = pp;
PhiPhi1 = PhiPhi;
I1 = I;
N1 = length(I);
clear nn pp PhiPhi I
load('1.mat')
N2 = length(I);
nn = [nn1,nn];
pp = [pp1,pp];
PhiPhi = [PhiPhi1,PhiPhi];
I = [I1,I];
for i=1:N1+N2
    hFactor_NO_DOP(i,1) = Vcp*pp(N/2,i)/Vce/nn(N/2,i);
    [q_NO(i,1),R(i),fn2(i),fn1(i),fp(i),P1_NO_DOP(i,1),P2_NO_DOP(i,1)] = F_QD_A2(Vce*nn(N/2,i)*c.ni,Vcp*pp(N/2,i)*c.ni);
    I_NO_DOP(i,1) = I(i);
end
clear nn pp PhiPhi I
load('22.mat')
nn1 = nn;
pp1 = pp;
PhiPhi1 = PhiPhi;
I1 = I;
N1 = length(I);
clear nn pp PhiPhi I
load('2.mat')
N2 = length(I);
nn = [nn1,nn];
pp = [pp1,pp];
PhiPhi = [PhiPhi1,PhiPhi];
I = [I1,I];
for i=1:N1+N2
    hFactor_P_DOP(i,1) = Vcp*pp(N/2,i)/Vce/nn(N/2,i);
    [q_P(i,1),R(i),fn2(i),fn1(i),fp(i),P1_P_DOP(i,1),P2_P_DOP(i,1)] = F_QD_A2(Vce*nn(N/2,i)*c.ni,Vcp*pp(N/2,i)*c.ni);
    I_P_DOP(i,1) = I(i);
end
clear nn pp PhiPhi I
load('33.mat')
nn1 = nn;
pp1 = pp;
PhiPhi1 = PhiPhi;
I1 = I;
N1 = length(I);
clear nn pp PhiPhi I
load('3.mat')
N2 = length(I);
nn = [nn1,nn];
pp = [pp1,pp];
PhiPhi = [PhiPhi1,PhiPhi];
I = [I1,I];
for i=1:12
    hFactor_N_DOP(i,1) = Vcp*pp(N/2,i)/Vce/nn(N/2,i);
    [q_N(i,1),R(i),fn2(i),fn1(i),fp(i),P1_N_DOP(i,1),P2_N_DOP(i,1)] = F_QD_A2(Vce*nn(N/2,i)*c.ni,Vcp*pp(N/2,i)*c.ni);
    I_N_DOP(i,1) = I(i);
end
figure(1)
plot(I_NO_DOP,P1_NO_DOP,I_P_DOP,P1_P_DOP,I_N_DOP,P1_N_DOP)
legend('NO-DOP','P-DOP','N-DOP')
figure(2)
plot(I_NO_DOP,P2_NO_DOP,I_P_DOP,P2_P_DOP,I_N_DOP,P2_N_DOP)
legend('NO-DOP','P-DOP','N-DOP')
figure(3)
plot(I_NO_DOP,P2_NO_DOP,I_P_DOP,P2_P_DOP,I_N_DOP,P2_N_DOP,I_NO_DOP,P1_NO_DOP,I_P_DOP,P1_P_DOP,I_N_DOP,P1_N_DOP)
legend('NO-DOP','P-DOP','N-DOP')
figure(4)
plot(I_NO_DOP,hFactor_NO_DOP,I_P_DOP,hFactor_P_DOP,I_N_DOP,hFactor_N_DOP)
legend('NO-DOP','P-DOP','N-DOP')