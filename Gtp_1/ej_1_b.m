clc; clear all;

% Ejercicio 1.b.

  L1 = 0;
  L2 = 2;
  N = 4;

  T_1 = 100;
  % T_2 := Neumann = q = 0
  q = 0;

  k = 1;
  c = 1;
  G = 0;
  dx = (L2 - L1)/(N-1);

  % Armar matriz K
  K = diag(-ones(N-1,1),-1) + diag((2 + c*dx^2/k)*ones(N,1)) + diag(-ones(N-1,1),1);

  % Armar vector f
  f = zeros(N,1);

  % Reemplazar en nodos interiores:
  f(2:N-1) = G/k*dx^2;

  % Aplicar condiciones de borde:
  f(1) = T_1; % Borde izq. : Dirichlet
  f(end) = G/k*dx^2 - 2*dx*q/k; % Borde der. : Neumann

  K(1, 1:2) = [1 0];
  K(end, end-1) = -2;

##  disp("Matriz K:")
##  K

##  disp("\nVector f:")
##  f

  disp("\nSolución por Diferencias Finitas:")
  T = K\f;
  T(end-2:end)

  % Analitico:
  T_a = @(x) (100 .* exp(-x).*(exp(2 .*x) + exp(4)))./(1 + exp(4));

  x_a = L1:dx:L2;
  disp(length(x_a))
  disp(x_a(2) - x_a(1))
  disp(dx)

  disp("\nSolución Analítica:")
  sol_a = T_a(x_a');
  sol_a(end-2:end)

  figure(1, 'name', "Grafico comparación")
  plot(x_a, T, 'linewidth', 1.5, 'r-')
  grid on; hold on;
  plot(x_a, sol_a, 'ko-', 'markersize', 5, 'markerfacecolor', 'b', 'linewidth', 1.5, 'color', 'c')

  ylabel("Temperatura")
  xlabel(" xi ")
  legend("MDF", "Analitico")
