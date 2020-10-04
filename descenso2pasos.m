function [W, H] = descenso2pasos(X,k)
% Metodo del descenso en dos pasos para encontrar una factorizacion no
% negativa de una matriz X
% Min norm(X-WH,'fro')^2
% s.a. W>=0, H>=0
%Llamado: [W, H] = descenso2pasos(X,k)
% In
% X.- matriz rxp no negativa.
% k.- es un numero natural tal que k<<r y k<<p
%
%Out
% W.- matriz rxk no negativa
% H.- matriz de kxp no negativa
%--------------------------------------------------------------------------
% Andres Cruz y Vera 155899
% Alexis Ayala Redon 156916
% Javier Montiel Gonzalez 159216
%--------------------------------------------------------------------------
[r, p] = size(X);
W = ones(r, k);
H = ones(k, p);
maxiter=4;

for i = 1:maxiter
    
    % Sea W fija
    % Resolvemos el problema 
    % Min norm(X-W*H,'fro')^2
    % s.a. H>=0
    % Minimizamos para la n-esima columna de X(*,n) y H(*,n)
    % Se resuelven p problemas cuadraticos
    for n = 1:p
        col_n_x = X(:, n);
        Q = W'*W; % dim k*k
        c = (-col_n_x'*W)';%  dim (1*r*r*k)' = k*1
        b = zeros(k,1);
        A = eye(k);
        [x, y, mu] = punintpc(Q,A,c,b);
        H(:,n) = x;   
    end
    
    % Sea H fija
    % Resolvemos el problema 
    % Min norm(X-W*H,'fro')^2
    % s.a. W>=0
    % Minimizamos para el m-esimo renglon de X(m,*) y W(m,*)
    % Se resuelven r problemas cuadraticos
    for m = 1:r
        renglon_m_x = X(m, :);%dim 1*p
        Q = H*H';%dim k*k
        c = (-renglon_m_x*H')'; %dim (1*p*p*k)' = k*1
        b = zeros(k,1);
        A = eye(k);
        [x, y, mu] = punintpc(Q,A,c,b);
        W(m,:) = x';   
    end
end

end

