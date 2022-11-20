function Y = addNoise(I,noiseLevel)
%ADDNOISE Summary of this function goes here
%   Detailed explanation goes here


Y=imnoise(I,'gaussian',0.0001,noiseLevel);
end

