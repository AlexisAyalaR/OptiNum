function [W, H] = descenso2pasos(X, k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[r, p] = size(X);
W = ones(r, k);
H = ones(k, p);
maxiter = 50;

for i = 1:maxiter
    %Resolvemos problema cuadrárico para W fija
    for n = 1:p
        col_n_x = X(:, n);
        Q = W'*W; % dim k*k
        c = (-col_n_x'*W)';
        b = zeros(k,1);
        A = eye(k);
        [x, y, mu] = punintpc(Q,A,c,b);
        H(:,n) = x;
        
    end
    %Resolvemos para H fija
    for m = 1:r
        renglon_m_x = X(m, :);%dim 1*p
        Q = H*H';%dim k*k
        c = (-renglon_m_x*H')'; %dim (1*p*p*k)' = k*1
        b = zeros(k,1);
        A = eye(k);
        [x, y, mu] = punintpc(Q,A,c,b);
        W(m,:) = x';
        
    end
    i
end

end

