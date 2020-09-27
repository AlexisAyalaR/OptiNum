function [x, y, mu] = punintpc(Q,A,c,b)
% Metodo de punto interior para el problema cuadr�tico
% Min   (0.5)* x' * Q * x + c�* x
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
% lambda.- vector en R^m con la aproximaci�n al multiplicador
%          de Lagrange asociado a la restricci�n de igualdad.
% mu.- vector en R^p con la aproximaci�n al multiplicador
%          de Lagrange asociado a la restricci�n de desigualdad.
% y.- vector en R^p con la variable de holgura en la restricci�n
%     de desigualdad.
%---------------------------------------------------------------------------
% Andr�s Cruz y Vera 155899
% Alexis Ayala Redon 156916
% Javier Montiel Gonzalez 159216
%--------------------------------------------------------------------------------
% Parametros iniciales
tol = 1e-08;       % Tolerancia a las condiciones necesarias de 1er orden
maxiter = 250;     % m�ximo n�mero de iteraciones permitidas
iter = 0;          % contador de las iteraciones
%-----------------------------------------------------------
n = length(c);     % dimensi�n de la variable principal
m = length(b);     % n�mero de restricciones 
%----------------------------------------------------------
% variables iniciales
mu = ones(m,1);
y = ones(m,1);

%Para obtener x0 resolvemos el situiente problema
epsilon = 1e-03;
e = ones(m,1);
Aeq = [A -1*eye(m)];
beq = b - epsilon*e;
x = linprog([e; zeros(n,1)],A,b,Aeq,beq);

%-----------------------------------------------
% vectores para graficaci�n
cnpo=[]; comp =[];
% Norma de las condiciones necesarias de primer orden
H =[Q*x + A'*lambda - F'*mu+c; A*x - b;F*x - z - d; mu.*z];
norma = norm(H);
disp('Iter      CNPO             tau ')
disp('-----------------------------------------')
while(norma > tol & iter < maxiter)
  % Resuelve el sistema lineal de Newton para la trayectoria central
    D = diag(mu./z);
    G = Q + F'*D*F;
    w = zeros(p,1);
    for k = 1:p
        w(k) = F(k,:)*x - d(k)-(tau/mu(k));
    end
    dg = Q*x + A'*lambda  -F'*mu+ c + F'*D*w;
    
    % sistema lineal
    K = [G  A'; A  zeros(m)];
    ld = -[dg ; A*x-b];
    y = K \ ld;
    %-------------------------------------------------------
    % Se calculan los pasos
    Dx = y(1:n);
    Dlambda = y(n+1:n+m);
    
    Dmu =-(D)*(F*Dx+w);
    Dz = -( (1./mu).*(z.*Dmu - tau) + z );
    %Dz =-inv(MatU)*(MatZ*Dmu+MatU*MatZ*vt-tau*vt);  
   %---------------------------------------------------------- 
    % Acorta el paso
    bt = []; gm = [];
    for k =1:p
        if (Dmu(k) < 0)
            bt = [bt; -(mu(k)/Dmu(k))];
        end
        if(Dz(k) < 0)
            gm = [gm; -(z(k)/Dz(k))];
        end
    end
    
    alfa = min([bt ; gm]);
    alfa =(0.9995)*min([1 alfa]);  
    %-----------------------------------------------------------
     % Nuevo punto
       x      = x + alfa*Dx;
       lambda = lambda + alfa*Dlambda;
       mu     = mu + alfa*Dmu;
       z      = z + alfa*Dz;
     %-------------------------------------------------------  
     % Nueva tau
        tau = (0.5)*(mu'*z)/p;
     %-------------------------------------------------------  
       %Condiciones necesarias de primer orden
       H=[Q*x+A'*lambda-F'*mu+c;A*x-b; F*x-z-d; mu.*z ];
       norma = norm(H);
       iter = iter + 1;
       cnpo =[cnpo norma];
       comp = [comp 2*tau];
       disp(sprintf('%3.0f  %2.8f  %2.8f',iter,norma,2*tau))
end

   semilogy([1:iter],cnpo,'r',[1:iter],comp,'b')
   title('Convergencia de puntos interiores')
   legend('CNPO', 'Complementaridad')
end
        
