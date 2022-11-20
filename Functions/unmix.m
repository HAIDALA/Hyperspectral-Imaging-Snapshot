%% 
function [G,F,iter,elapse,HIS ] = unmix(I,r,Ginit,Finit,alpha)
% This function is to do unmxing with Nesterov NMF solver
% alpah parameter is to control sum-to-one constraint

% Created by Kinan ABBAS on January 2022
% Last update date Oct 5 2022
%<Outputs>
%        G : Obtained basis matrix (m x r).
%        F : Obtained coefficients matrix (r x n).
%        iter : Number of iterations.
%        elapse : CPU time in seconds.
%        HIS : (debugging purpose) History of computation,
%               niter - total iteration number spent for Nesterov's optimal
%               gradient method,
%               cpus - CPU seconds at iteration rounds,
%               objf - objective function values at iteration rounds,
%               prjg - projected gradient norm at iteration rounds.

[n1,n2,n3]=size(I);

%Creating spectral measurment matrix
% Y size is (n1*n2)*n3. Each column is a single image under different
% wavelenght
Y=zeros(n1*n2,n3);
for i=1:n3
    temp=I(:,:,i);
    Y(:,i)=reshape(temp,[n1*n2,1]);
end

%Sum to one constraint
Y1=[Y,alpha*ones(n1*n2,1)]; Finit1=[Finit,alpha*ones(r,1)];
% Y1=Y;
% Finit1=Finit;
% [ F , G , iter,elapse,HIS]=NeNMF(Y',r,'TYPE','L1R','BETA',0,'MAX_ITER',1000,'MIN_ITER',10,'W_INIT',Finit','H_INIT',Ginit');
%  [ F , G , iter,elapse,HIS]=NeNMF(Y',r,'TYPE','L1R','BETA',0,'MAX_ITER',100000,'MIN_ITER',10000);
% [ F , G , iter,elapse,HIS]=NeNMF(Y',r,'SCALLING',false,'MAX_ITER',1000,'MIN_ITER',10,'W_INIT',Finit','H_INIT',Ginit');
[ F , G , iter,elapse,HIS]=NeNMF(Y1',r,'SCALLING',false,'MAX_ITER',1000,'MIN_ITER',10,'W_INIT',Finit1','H_INIT',Ginit');

F=F(1:n3,:)';
G=G';

