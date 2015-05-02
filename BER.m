function [ y ] = BER( x )
%BER Summary of this function goes here
%   Detailed explanation goes here
    if (x ~= 0)
        y = x/(exp(x)-1);
    else
        y = 1;
    end
end

