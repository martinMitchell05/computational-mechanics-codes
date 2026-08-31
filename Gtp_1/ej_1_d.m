% Ejercicio 1. d
clc; clear all;

  L1 = 0;
  L2 = 1;

  T_0 = 10;
  % T_3 = robin ---> aplico aproximación O(dx^2) con nodo ficticio
  h = 0.2;
  T_inf = 50;

  N = 4; % cantidad de puntos
  dx = (L2 - L1)/(N-1);

  k = 1;
  G = 50;
  c = 1;

  % Armar matriz 'K'
  K = diag(-ones(N-1,1), -1) + diag((2 + (c*dx^2)/k)*ones(N,1)) + diag(-ones(N-1,1), 1)

  % Armar vector 'f'
  f = zeros(N,1);
  f(2:N-1) = G/k * (dx)^2;

  % Aplicar condición de borde
  f(1) = T_0;
  f(end) = G/k * (dx)^2 + 2 * dx * h * T_inf/k;
  disp("\n")
  f

  % Aplicar condición de borde
  K(1,1:2) = [1 0];
  K(end, end-1:end) = [-2, (2*(1 + dx*h) + c*dx^2)];
  disp("\nMatriz 'K' final:")
  disp(K)

  % Resolver sistema:
  T = K\f;
  disp("\nSolución a la ecuación diferencial por Diferencias Finitas:")
  disp(T)


  % Solución analítica:
  x_a = linspace(L1, L2, N)';
  T_a = @(x) -36.6897*exp(-x) - 3.3103*exp(x) + 50;

  disp("\nSolución analítica:")
  disp(T_a(x_a))

  figure(1, 'name', "Grafico comparación")
  plot(x_a, T, 'linewidth', 1.5, 'r-')
  grid on; hold on;
  plot(x_a, T_a(x_a), 'ko-', 'markersize', 5, 'markerfacecolor', 'b', 'linewidth', 1.5, 'color', 'c')

  ylabel("Temperatura")
  xlabel(" xi ")
  legend("MDF", "Analitico")


