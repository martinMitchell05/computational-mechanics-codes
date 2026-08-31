clc; clear all;
% Ejercicio 1.c

  L1 = 1;
  L2 = 5;

  % T_1 ---> Neumann = 2
  T_2 = 0;
  k = 1;
  c = 0;
  G = @(x) 100 .* (x - 3).^2;

  N = 12;
  dx = (L2 - L1)/(N-1);
  xi = (L1:dx:L2);

  % Armar matriz 'K'
  K = diag(-ones(N-1,1), -1) + diag((2 + c/k*(dx^2))*ones(N,1)) + diag(-ones(N-1,1), 1)

  % Armar vector 'f'
  f = zeros(N,1);
  f(2:N-1) = G(xi(2:N-1))./k.*(dx^2)

  % Aplicar condiciones de borde
  K(end, end-1:end) = [0 1];
  % f(end) es 0

  K(1, 1:2) = [(2 + c/k*(dx^2)), -2];
  f(1) = (dx*(G(xi(1))*dx - 4))/k;

  disp("\nMatriz final K:")
  K

  disp("\nVector final f:")
  f

  disp("\nResolución por Diferencias Finitas:")
  T = K\f

  % Solución analítica:
  a = -1350;
  b = 1906;
  c = 2345;

  T_a = @(x) (-25.*x.^4 + 300.*x.^3 + a.*x.^2 + b.*x + c)./3;

  disp("\nSolución analítica:")
  disp(T_a(xi'))

  figure(1, 'name', "Grafica comparación")
  plot(xi, T, 'linewidth', 1.5, 'color', 'r')
  grid on; hold on;
  plot(xi, T_a(xi), 'ko-', 'markersize', 5, 'linewidth', 1.5, 'color', 'b')

  xlabel("  xi  ")
  ylabel("Temperatura")
  legend("MDF", "Analitico")





