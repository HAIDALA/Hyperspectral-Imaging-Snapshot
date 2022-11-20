function [F] = Scale(F,HighestNumber)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
maxi=max(max(F));
F=F/maxi;
F=F*(HighestNumber);
% F=round(F);
F=F/10;
end