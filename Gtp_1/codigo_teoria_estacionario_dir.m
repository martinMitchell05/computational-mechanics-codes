    % EcPoissonDirichletFD.m


clear all; clc;
% programa para resolver la ecuacion de Poisson
% con condiciones Dirichlet
% Simil ecuacion de calor con Fuente

% Metodo: Diferencias Finitas

% -kphi''(x) = Q(x) para todo x, tal que 0<x<L
% phi(0) = a; phi(L) = b;

% parametros del problema

L=1;  k=1; a=2; b=-3;
% para diferentes fuentes:
 Q = @(x) 100*exp(-20*(x-L/3).^2);
% Q = @(x) 300*(0.5<x<0.7);
% Q = @(x) 100*(sin(5*pi*x));
% Q = @(x) 300*(0.3<x<0.4)-(-100*(0.7<x<0.8));
% Q = @(x) 300*(0.3<x<0.4)+300*(sin(7*pi*x));
%Parametros de metodo de discretizacion y resolucion

N=30;

% Armado de la matriz
% paso
h = L/N;
uno=ones(N+1,1);
diagonales =[-1*uno 2*uno -1*uno];
matriz=spdiags(diagonales, [-1 0 1], N+1, N+1);
% Arreglo primera fila de acuerdo condiciones de contorno
matriz(1,[1:2]) = [1 0];
% Arreglo ultima fila de acuerdo condiciones de contorno                                                                                              +
matriz(N+1, [N:N+1]) = [0 1];

% Armamos el lado derecho
x=linspace(0,L,N+1)';
F = [a; (h^2/k).*Q(x(2:N)); b];

% Resolvemos:
phi = matriz\F; % método directo de Octave

% ANALIZAMOS LA MATRIZ
%=====================================
rhoMatriz = max(abs(eig(matriz)))
##edd(matriz)
##analisisMatriz(matriz,F)
%=====================================
% graficamos
figure(1)
plot(x,phi,'b*-')
grid on; grid minor
title('Ecuacion de Poisson - Condiciones Dirichlet')
xlabel('longitud dominio: X')
ylabel('Solucion: Temperatura')
