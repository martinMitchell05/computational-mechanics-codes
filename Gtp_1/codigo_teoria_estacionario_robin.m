% EcPoissonDirichletFD.m
% programa para resolver la ecuacion de Poisson
% con condiciones Dirichlet y Mixtas (robin) o Neumann
% Simil ecuacion de calor con Fuente
clear all; clc;
% Metodo: Diferencias Finitas

% -kphi''(x) = Q(x) para todo x, tal que 0<x<L
% phi(0) = a; k*phi'(L) + H1*phi(L)=H2*phi_E

% parametros del problema

L=1;  k=1; a=2; H1=1; H2=1; phi_E=0;
% para diferentes fuentes:
Q =@(x) 100*exp(-20*(x-0.3).^2);
% Q = @(x) 300*(x>0.7);
% Q = @(x) 300*(0.5<x<0.7);
% Q = @(x) 100*(sin(5*pi*x));
% Q = @(x) 300*(0.3<x<0.4)-(-100*(0.7<x<0.8));
% Q = @(x) 300*(0.3<x<0.4)+300*(sin(7*pi*x));
%Parametros de metodo de discretizacion y resolucion

N=30;

% Armado de la matriz
% paso
h = L/N;
uno=ones(N+2,1);
diagonales =[-1*uno 2*uno -1*uno];
matriz=spdiags(diagonales, [-1 0 1], N+2, N+2);
% Arreglo primera fila de acuerdo condiciones de contorno
matriz(1,[1:2]) = [1 0];
% Arreglo ulrima fila de acuerdo condiciones de contorno                                                                                              +
matriz(N+2, N:N+2) = [-1 2*h*H1/k 1];

% Armamos el lado derecho
x=linspace(0,L,N+1)';
F = [a; (h^2/k).*Q(x(2:N+1)); 2*h*H2*phi_E/k];

% Resolvemos:
phi = matriz\F;

% Elimino el punto ficticio
phi(N+2) = []; % otra manera: phi =phi(1:N+1);
% graficamos
figure(2)
plot(x,phi,'b*-')
grid on; grid minor
title('Ecuacion de Poisson - Condiciones Mixtas')
xlabel('longitud dominio: X')
ylabel('Solucion: Temperatura')
