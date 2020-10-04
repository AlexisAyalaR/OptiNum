clear all; clc;
%--------------------------------------------------------------------------
% Andres Cruz y Vera 155899
% Alexis Ayala Redon 156916
% Javier Montiel Gonzalez 159216
%--------------------------------------------------------------------------
%Cargamos la imagen del payaso de MATLAB
load clown;
%La transformamos a una escala de grises
colormap('gray');
%Desplegamos la imagen original X
figure (1)
image(X);

%Buscamos factorizacion no negativa de X
[W, H] = descenso2pasos(X, 60);
%Calculamos el error de aproximación
norm(X-W*H,'fro')
%Desplegamos la imagen obtenida
image(W*H);