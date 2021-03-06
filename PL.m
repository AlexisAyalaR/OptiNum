%Para obtener x0 resolvemos el situiente problema
%Para obtener x0 resolvemos el situiente problema
clear all; clc;
Q=[1 0 0; 0 1 0; 0 0 1];
c=[-1;-1;-1];
A=[1 0 0; 0 1 0; 0 0 1;-1 0 0; 0 -1 0; 0 0 -1];
b=[0; 0; 0; -1; -1; -1];
[x1, y, mu] = punintpc(Q,A,c,b);

n=length(c);
m=length(b);
epsilon = 1e-09;
e = ones(m,1);
Aeq = [-eye(m) A];
beq = b - epsilon*e;
x = linprog([e; zeros(n,1)],-[eye(m,m+n);zeros(n,m+n)],-zeros(m+n,1),Aeq,beq);