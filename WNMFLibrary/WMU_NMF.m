function [ G , F ] = WMU_NMF( W , X , G , F , Iter_max )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %
%   Weighted extension (from Guillamet) of Lee & Seung Multiplicative update rules %
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %
%  solve:    argmin{G,F} || W.*(X - G.F) ||^2_Fro                %
%                                                                %
%            s.t. G , F >= 0                                     %
%                                                                %
%  input: W = m*n matrix of weights                              %
%         X = m*n matrix of data                                 %
%         Ginit = m*r matrix (initialization for G)              %
%         Finit = r*n matrix (initialization for F)              %
%                                                                %
%  output: estimated G and F                                     %
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   author : Clément Dorffer                                     %
%   date : 15/05/2016                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initf = norm( W.*(X-G*F) , 'fro' )^2;
W2 = W.^2;
WX = W2.*X;
tic
for i = 1 : Iter_max    
    % updating G
    G = Updt_G( W2 , WX , G , F );    
    % updating F
    F = Updt_F( W2 , WX , G , F );    
end
t = toc;
f = norm( W.*(X-G*F) , 'fro' )^2;
fprintf('\n###   Elapse time: %d sec.\n###   Initial objective value: %d\n###   Objective value: %d \n' , t , initf , f );
end



function [ F ] = Updt_F( W , WX , G , F )
%%%%%%%%%%%%%%%%%%%%%
% update rule for F %
%%%%%%%%%%%%%%%%%%%%%
F = F.*((G'*WX)./secu_plus(G'*(W.*(G*F)),eps));
F(isnan(F)) = 0;
end



function [ G ] = Updt_G( W , WX , G , F )
%%%%%%%%%%%%%%%%%%%%%
% update rule for G %
%%%%%%%%%%%%%%%%%%%%%
G = G.*((WX*F')./secu_plus((W.*(G*F))*F',eps));
G(isnan(G)) = 0;
end



function [toto] = secu_plus(tutu,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Matthieu Puigt
% Date: 20/04/2015
% Goal: Security in the NMF procedure which project the negative data to
% epsilon (small user-defined threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin<2)
    s=1.e-8;
end
tutu(isnan(tutu)) = 0;
toto=max(tutu,s);
end
