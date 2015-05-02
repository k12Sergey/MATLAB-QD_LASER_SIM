N = 350;
[ qd ] = QD_parameters();
[ c ] = constants();
[ Nc,Nv,ni1,eps,P,Eg,taun,taup,mun,mup,Xi,C(:,1),Theta_n1,Theta_p1,X,h,b,L,D1 ] = StructureProfile( N );
Nqd = fix(qd.l/h(1));   % width of QD layer  fix(qd.l/h(1))
Nqd2 = N/2-fix(Nqd/2);  % place for QD layer, in a middle
Vce = qd.gen/Nc(Nqd2)/fermi(1/2,-(qd.dEc1-qd.en2)/c.Vt);
Vcp = qd.gep/Nv(Nqd2)/fermi(1/2,-(qd.dEv1-qd.ep)/c.Vt);

load('new_111111.mat');
N1 = length(I);
for i=1:N1
     hFactor_NO(i,1) = Vcp*pp(N/2,i)/Vce/nn(N/2,i);
    [q_NO(i,1),R(i),fn2(i),fn1(i),fp(i),P1_NO_DOP(i,1),P2_NO_DOP(i,1)] = F_QD_A2(Vce*nn(N/2,i)*c.ni,Vcp*pp(N/2,i)*c.ni);
end