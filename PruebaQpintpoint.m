%Prueba del método punintpc
%Min (1/2)*[(x1-2)^2+(x2-2)^2+(x3-2)^2)
%s.a. 0<=xi<=1
%El óptimo se alcanza en (1,1,1)'
Q=[1 0 0; 0 1 0; 0 0 1];
c=[-1;-1;-1];
A=[1 0 0; 0 1 0; 0 0 1;-1 0 0; 0 -1 0; 0 0 -1];
b=[0; 0; 0; -1; -1; -1];
[x, y, mu] = punintpc(Q,A,c,b);