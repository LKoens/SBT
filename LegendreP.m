function [y] = LegendreP(n,x)
%LEGENDREP Summary of this function goes here
%   Detailed explanation goes here

L1= legendre(n,x);
if n==0
    y=L1;
else
y =reshape(L1(1,:,:),size(x));
end
end

