function [ F,Rad_R,Rad_L ] = RP( N,l1,l2,l3,n,p,X,V,t )
%   Right side of system with a.u.
%   If t = 1 HOLES, t = 2 ELECTRONS
    [ c ] = constants();
    b1 = l1.b;
    b2 = l1.b+l2.b;
    if(t == 1)                      % for holes
        Norm = muh(2)*c.Vt/c.Ld^2;        
    elseif(t == 2)                  % for electrons
        Norm = mun(2)*c.Vt/c.Ld^2;        
    end
    Nm = 1/Norm;
    ni1 = l1.ni/c.ni*exp(-V/c.Vt);
    ni2 = l2.ni/c.ni*exp(-V/c.Vt);
    ni3 = l3.ni/c.ni*exp(-V/c.Vt);
    for i = 1:N-1
        if(X(i) < b1)
            F(i) = Nm*1/(l1.tau_p*(n(i)+ni1)+l1.tau_n*(p(i)+ni1));
            Rad_R(i) = Nm*c.ni*l1.B*ni1^2;
            Rad_L(i) = Nm*c.ni*l1.B;
        elseif(X(i) <= b2)
            F(i) = Nm*1/(l2.tau_p*(n(i)+ni2)+l2.tau_n*(p(i)+ni2));
            Rad_R(i) = Nm*c.ni*l2.B*ni2^2;
            Rad_L(i) = Nm*c.ni*l2.B;
        else
            F(i) = Nm*1/(l3.tau_p*(n(i)+ni3)+l3.tau_n*(p(i)+ni3));
            Rad_R(i) = Nm*c.ni*l3.B*ni3^2;
            Rad_L(i) = Nm*c.ni*l3.B;
        end
    end
%     F= -F;
    F = zeros(N,1);
    Rad_R = zeros(N,1);
    Rad_L = zeros(N,1);
end

