function [ y ] = Destribution( Phi,ni,t,k,Phi_n,Phi_p,Theta_n,Theta_p,Vr )
% DIstribution function statistics
%   Ef = 0 Bolzman statistics, Phi = Phi/c.Vt you need to provide
    [ c ] = constants();
    switch t
        case 1  % n
            switch k
                case 1  % function
                    y = ni.*exp(Theta_n+(Phi-Vr)/c.Vt-Phi_n);
                case 2  % derivative Phi
                    y = ni.*exp(Theta_n+(Phi-Vr)/c.Vt-Phi_n);                
            end
         case 2  % p
            switch k
                case 1  % function
                    y = ni.*exp(Theta_p+(-Phi+Vr)/c.Vt+Phi_p);
                case 2  % derivative Phi
                    y = -ni.*exp(Theta_p+(-Phi+Vr)/c.Vt+Phi_p);
            end
    end

end

