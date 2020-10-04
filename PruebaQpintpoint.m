clear all; clc;
%Prueba 1 del método punintpc
%Min (1/2)*[(x1-2)^2+(x2-2)^2+(x3-2)^2)
%s.a. 0<=xi<=1
%El óptimo se alcanza en (1,1,1)'
Q=[1 0 0; 0 1 0; 0 0 1];
c=[-1;-1;-1];
A=[1 0 0; 0 1 0; 0 0 1;-1 0 0; 0 -1 0; 0 0 -1];
b=[0; 0; 0; -1; -1; -1];
[x, y, mu] = punintpc(Q,A,c,b);

%Prueba 2 del método puntintpc
% Min (1/2)*[2*x1^2+2*x2^2]
% s.a. xi>=0
Q1=[2 0;0 2];
c1=[0 ;0];
A1=eye(2);
b1=[0;0];
[x1, y1, mu1] = punintpc(Q1,A1,c1,b1);
%x = quadprog(Q1,c1,-A1,-b1);
