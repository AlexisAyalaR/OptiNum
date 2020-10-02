function [W, H] = descenso2pasos(X, k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[r, p] = size(X);
W = ones(r, k);
H = ones(k, p);
maxiter = 300;

for i = 1:maxiter
    for n = 1:p
        col_n_x = X(:, n);
        Q = W'*W;
        c = -col_n_x'*W;
        b = zeros(r,1);
        [x, y, mu] = punintpc(Q,A,c,b);
        H(:,n) = x;
    end

end

end

