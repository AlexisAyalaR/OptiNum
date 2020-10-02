function [x, y, mu] = punintpc(Q,A,c,b)
% Metodo de punto interior para el problema cuadr�tico
% Min   (0.5)* x' * Q * x + c'* x
% s.a.   A * x >= b
%
%  Llamado: function [x, y, mu] = punintpc(Q,A,c,b)
% In
% Q.- matriz nxn sim�trica y positiva definida
% A.- matriz mxn con m <= n y rango(A) = m..
% b.- vector columna en R^m .
% c.- vector columna en R^p .
%
%Out
% x.- vector en R^n con la aproximaci�n del m�nimo local.
% lambda.- vector en R^m con la aproximaci�n al multiplicador de Lagrange 
% asociado a la restricci�n de igualdad.
% mu.- vector en R^p con la aproximaci�n al multiplicador de Lagrange 
% asociado a la restricci�n de desigualdad.
% y.- vector en R^p con la variable de holgura en la restricci�n de 
% desigualdad.
%--------------------------------------------------------------------------
% Andr�s Cruz y Vera 155899
% Alexis Ayala Redon 156916
% Javier Montiel Gonzalez 159216
%--------------------------------------------------------------------------
% Parametros iniciales
tol = 1e-08;       % Tolerancia a las condiciones necesarias de 1er orden
maxiter = 250;     % m�ximo n�mero de iteraciones permitidas
iter = 0;          % contador de las iteraciones
%--------------------------------------------------------------------------
n = length(c);     % dimensi�n de la variable principal
m = length(b);     % n�mero de restricciones 
%--------------------------------------------------------------------------
% variables iniciales
x = zeros(n,1);
mu = ones(m,1);
y = ones(m,1);
e = ones(m,1);
nu = (0.5)*(mu'*y)/m;

% Norma de las condiciones necesarias de primer orden
H =[Q*x - A'*mu+c; A*x- y - b; mu.*y];
norma = norm(H);

while(norma > tol && iter < maxiter)
    Y=diag(y);
    U=diag(mu);
    D=diag(mu./y); %Y inversa por U
    
    % Matriz Jaconiaba del Lagrangiano
    rx=Q*x-A'*mu+c; %nx1
    ry=A*x-y-b;     %mx1
    rmu=Y*U*e-nu*e; %mx1
    
    % Resolvemos el sistema 
    B=Q+A'*D*A;
    d=-(rx+A'*D*ry+A'*(rmu./y));
    Dx= d\B;
    
    % Obtenemos los valores de Dy y Dmu
    Dy=A*Dx'+ry;
    Dmu=-(D*Dy-rmu./y);
    
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
       x      = x + alfa*Dx';
       mu     = mu + alfa*Dmu;
       y      = y + alfa*Dy;
     %---------------------------------------------------------------------  
     % Nueva nu
       nu = (0.5)*(mu'*y)/m;
     %---------------------------------------------------------------------
     %Condiciones necesarias de primer orden
       H =[Q*x - A'*mu+c; A*x - y - b; mu.*y];
       norma = norm(H);
       iter = iter + 1;
end
end