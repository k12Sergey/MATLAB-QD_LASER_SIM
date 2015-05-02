clc;
[ c ] = constants();
[ qd ] = QD_parameters();
h = 0.11;
gcn = linspace(1e9,1e12,100);
gamma1 = 0.6;
gamma2 = gamma1/1.8;
for i=1:100
    [ q,R,fe2,fe1,fh,P1(i),P2(i),J(i) ] = F_QD_A22(gcn(i),h*gcn(i),gamma1,gamma2);
    [gcn(i),P1(i),P2(i)];
    I(i) = R*qd.l*c.q*qd.Ns;
end
plot(J,P1,J,P2)
