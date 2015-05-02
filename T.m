%% TEST
for j = 1:Nn2-1
    for i = 2:N
        E(i,j) = -(Phi2(i+1,j)-Phi2(i-1,j))/2/h(i)/c.Ld;         
    end
    [Ec(j,:),n12(j,:)] = Potential(E(800,j),c.kb,c.T,c.q,c.eps,c.eps0,l1.N,l1.Nc);    
    [Ec2(:,j),Ev2(:,j)] = BandD( N,Phi2(:,j),X );
end
figure(1)
plot(X(1:N),E(:,3),X(1:N),E(:,4),X(1:N),E(:,5),X(1:N),E(:,6),X(1:N),E(:,7))
figure(2)
plot(X,Ec2(:,3),X,Ec2(:,4),X,Ec2(:,5),X,Ec2(:,6),X,Ec2(:,7),X,Ev2(:,3),X,Ev2(:,4),X,Ev2(:,5),X,Ev2(:,6),X,Ev2(:,7))