function [ Q,R,fe2,fe1,fh ] = F_QD_A_correct(gec,ghc)
[ qd ] = QD_parameters();
gee=qd.gen;
ghe=qd.gep;
ge21=qd.ge21;
ge12=qd.ge12;
gamma1=qd.gamma1;
gamma2=qd.gamma2;
ksi=ge12/ge21;

%NL option
fe2=gec/(gec+gee);
fh=ghc/(ghc+ghe);
fe1=fe2/(ksi+(1-ksi)*fe2);
P1=0;P2=0;
I1=fe1+fh-1;
I2=fe2+fh-1;
if((I1>gamma1) || (I2>gamma2))
	%NL option failed
	%try GS option
	B1=(gec-1.5*(ghc-(1+gamma1)*(ghc+ghe)))/(gec+gee);
	A1=1.5*(ghc+ghe)/(gec+gee);
	a=A1*(ge21-ge12);
	b=-(ge21*(A1+B1)+ge12*(1-B1)+3/2*(ghc+ghe));
	c=ge21*B1-3/2*ghc+3/2*(ghc+ghe)*(1+gamma1);
	d=sqrt(b^2-4*a*c);
	fe1=(-b-d)/2/a;
	fe2=B1-A1*fe1;
	fh=1+gamma1-fe1;
	P1=6*(ghc-(ghc+ghe)*fh);P2=0;
	I2=fe2+fh-1;
	if(I2>gamma2)
		%means ES or ES+GS
		fe2=(gamma2*ghc+ghe*(1+gamma2)+2/3*gec)/...
			(ghc+ghe+2/3*(gec+gee));
		fh=1+gamma2-fe2;
		fe1=fe2/(ksi+(1-ksi)*fe2);
		if((fe1+fh-1)>gamma1)
			%GS+ES
			fe1=1+gamma1-fh;
			P1=4*(ge21*fe2*(1-fe1)-ge12*fe1*(1-fe2));
		else
			P1=0;
		end
		P2=6*(ghc-(ghc+ghe)*fh)-P1;
	end
end

good1=2*(ge21*fe2*(1-fe1)-ge12*fe1*(1-fe2))-P1/2;
good2=gec-(gec+gee)*fe2-P1/4-P2/4;
good3=ghc-(ghc+ghe)*fh-(P1+P2)/6;
SOPLI=1e-5*ghe;
if(abs(good1)>SOPLI || abs(good2)>SOPLI || abs(good1)>SOPLI)
	panic_mode='ON'
end
R = 6*(ghc*(1-fh)-ghe*fh);
Q=6*fh-4*fe2-2*fe1;
end