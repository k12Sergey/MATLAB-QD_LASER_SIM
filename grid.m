function [ X,h ] = grid( b,N )
% Buid up the grid
%   For dimlicity use 
%     TMP = load('405.txt');
%     X = TMP(:,2);
%     X = sort(X)*1e-4;
%     N = length(X);
    X = linspace(0,b,N+1);
    h = X(2:N+1)-X(1:N);
% GRID
%     X = linspace(0,1,N+1)*b;
%     h = ones(N+1,1)*(X(3)-X(2));
end

