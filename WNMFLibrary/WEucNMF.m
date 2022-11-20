function [T, W, H, Func] = WEucNMF(V, W0, H0, MaxIter,Gamma)

%% [T, W, H, Func] = WEucNMF(V, W0, H0, MaxIter,Gamma)
%
% Authers:Jiao Wei, Can Tong, Bingxue Wu, Shouliang Qi, Yudong Yao, and Yueyang Teng
%Email: 2501554500@qq.com
% Entropy weighted non-negative matrix factorization (EWNMF) implementation
% solved by an iterative algorithm.
%
% Given a non-negative matrix V, a optimal hyperparamter Gamma, find non-negative matrix factors W and H 
% such that V approx. W*H, i.e. solving the following optimization problem:
% min_(T,W,H)=T.*||V-WH||^F_2 + Gamma*T.*ln(T)
%% inputs
%    V: data matrix (m x n)
%    W0, H0: initializations for W and H, W (m x k) is called base matrix , H (k x n) is
%    called representation matrix.
%    MaxIter: number of iterations
%    Gamma: a hyperparamper

%% outputs
%    W, H: factorization such that V \approx W*H
%    T: an optimizable weight matrix
%    Func: iterative error

% Test whether the dimensions of W0, H0 are correct
fprintf('WEucNMF\n')
if size(V) ~= size(W0*H0)
    fprintf('incorrect size of W0 or H0\n')
end

% Initialization of W, H, T, Func
[RowNum, ColNum] = size(V);
W = W0;
H = H0;
Func = zeros(MaxIter, 1);
T = ones(RowNum,ColNum) / RowNum;

% Iteration
for i = 1:MaxIter              
    
    % Update to W
    TV = T .* V;
    TWH = T .* (W * H);
    W = W .* ( TV * H') ./ ( TWH * H' +eps);
    
    %Update to H
    TWH = T .* (W * H);
    H = H .* ( W' * TV ) ./ (W' * TWH + eps);
    
    %Update to T
    T = exp( - (V - W * H).^2 / Gamma );
    T = T ./ repmat(sum(T) + eps, RowNum, 1);     
    
    %Iteration error
    Func(i) =  sum(sum(T.*(V - W * H).^2)) + Gamma * sum(sum(T.*log(T+eps)));
end

return;

