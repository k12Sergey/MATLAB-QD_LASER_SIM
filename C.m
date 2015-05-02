function [ C ] = C( N,X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [ l1,l2,l3,het,b ] = structure_parameters();
    [ c ] = constants();
    b1 = l1.b;
    b2 = l1.b+l2.b;
    for i=1:N+1
        if(X(i) <= b(1))
           C(i) = l1.N;         % Nd
        elseif(X(i) <= b(1)+b(2))
           C(i) = l2.N;
        elseif(X(i) <= b(1)+b(2)+b(3))
           C(i) = l2.N;
        elseif(X(i) <= b(1)+b(2)+b(3)+b(4))
           C(i) = l2.N;
        elseif(X(i) <= b(1)+b(2)+b(3)+b(4)+b(5))
           C(i) = l2.N;
        elseif(X(i) <= b(1)+b(2)+b(3)+b(4)+b(5)+b(6))
           C(i) = l2.N;
        elseif(X(i) <= b(1)+b(2)+b(3)+b(4)+b(5)+b(7))
           C(i) = -l3.N;        
        end
    end
    C = C;
end

