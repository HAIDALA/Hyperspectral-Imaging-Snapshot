function Y = addNoise(I,noiseLevel)
%ADDNOISE Summary of this function goes here
%   Detailed explanation goes here

[n1,n2]=size(I);
Y=I+noiseLevel*rand(n1,n2);
end

