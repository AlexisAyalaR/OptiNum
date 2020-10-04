function [x, y, mu] = punintpc(Q,A,c,b)
% Metodo de punto interior para el problema cuadratico
% Min   (0.5)* x' * Q * x + c'* x
% s.a.   A * x >= b
%
% Llamado: [x, y, mu] = punintpc(Q,A,c,b)
% In
% Q.- matriz nxn simetrica y positiva definida
% A.- matriz mxn con m <= n y rango(A) = m..
% b.- vector columna en R^m .
% c.- vector columna en R^p .
%
%Out
% x.- vector en R^n con la aproximacion del minimo local.
% mu.- vector en R^p con la aproximacion al multiplicador de Lagrange 
% asociado a la restriccion de desigualdad.
% y.- vector en R^p con la variable de holgura en la restriccion de 
% desigualdad.
%--------------------------------------------------------------------------
% Andres Cruz y Vera 155899
% Alexis Ayala Redon 156916
% Javier Montiel Gonzalez 159216
%--------------------------------------------------------------------------
% Parametros iniciales
tol = 1e-08;       % Tolerancia a las condiciones necesarias de 1er orden
maxiter = 250;     % maximo numero de iteraciones permitidas
iter = 0;          % contador de las iteraciones
%--------------------------------------------------------------------------
n = length(c);     % dimension de la variable principal
m = length(b);     % numero de restricciones 
%--------------------------------------------------------------------------
% variables iniciales
mu = ones(m,1);
y = ones(m,1);
e = ones(m,1);
Y= diag(y);
Yinv= diag(1./y);
U= diag(mu);

%Punto inicial
x= ones(n,1);

%Calculamos eta
eta = (0.5)*(mu'*y)/m;


% Norma de las condiciones necesarias de primer orden
H =[Q*x - A'*mu+c; A*x- y - b; Y*U*e];
norma = norm(H);

while(norma > tol && iter < maxiter)    
    % Condiciones perturbadas de KKT
    rx=Q*x-A'*mu+c; %nx1
    ry=A*x-y-b;     %mx1
    rmu=Y*U*e-eta*e; %mx1
    
    % Resolvemos el metodo de Newton 
    B=Q+A'*Yinv*U*A;
    d=-(rx+A'*Yinv*U*ry+A'*Yinv*rmu);
    Dx= B\d; 
    
    % Obtenemos los valores de Dy y Dmu
    Dy=A*Dx+ry;
    Dmu=-Yinv*(U*Dy+rmu);
    
    % Acortamos el paso
    bt = []; gm = [];
    for k =1:m
        if (Dy(k) < 0)
            gm = [gm; -(y(k)/Dy(k))];
        else
            gm = [gm; 1];
        end
        if(Dmu(k) < 0)
            bt = [bt; -(mu(k)/Dmu(k))];
        else
            bt= [bt; 1];
        end
    end
    
    alfa = min([bt ; gm]);
    alfa =(0.9995)*min([1 alfa]);  
    %----------------------------------------------------------------------
    % Nuevo punto
    x = x + alfa*Dx;
    mu = mu + alfa*Dmu;
    y = y + alfa*Dy;
    
    Y=diag(y);
    U=diag(mu);
    Yinv=diag(1./y); 
    %----------------------------------------------------------------------  
    %Nueva eta
    eta = (0.5)*(mu'*y)/m;
    %---------------------------------------------------------------------
    %Condiciones necesarias de primer orden
    H =[Q*x - A'*mu+c; A*x - y - b; Y*U*e];
    norma = norm(H);
    iter = iter + 1;
end
end