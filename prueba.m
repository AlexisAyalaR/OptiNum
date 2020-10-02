clear;
clc;
k = 3;
W_sol = [ 9 9 9; 7 7 7];
H_sol = [2 2; 1 1; 9 9];
X = W_sol*H_sol;
[W, H] = descenso2pasos(X,k,H_sol);
norm(X-W*H, 'fro')