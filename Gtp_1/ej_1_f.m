% Ejercicio 1.f
clc; clear all;

  L1 = 0;
  L2 = 1;
  T_1 = 0;

  N = 32;
  % T_2 : Robin
  h = 2;
  T_inf = 10;

  cp = 2;
  ro = 1;
  k = 2;
  c = 2;
  G = 75*ones(N,1);
  model = {k, c, ro, cp, G};

  dx = (L2 - L1)/(N-1);
  lambda = 0.9;
  dt = 0.05;
  tol = 1e-7;
  maxIt = 3000;

  xnode = L1:dx:L2;
  cb = [1 T_1 0; 3 h T_inf];

  % Resolucion estacionaria:
  et = [0 0 0 0 0];
##  et = [0 0 0 0];

  disp("Resolución por Diferencias Finitas:\n")
  tic
  T = difFinitas(xnode, model, cb, et);
  t_estacionaria = toc

  figure(1)
  plot(xnode, T, 'ko-', 'markersize', 7, 'linewidth', 1.5, 'color', 'b', 'DisplayName', 'Estacionaria')
  grid on; hold on;

  printf("Temperatura en nodo final x=L2 (Estacionario): T = %.5f\n", T(end));

  % Resolucion explicita:
  et = [1 maxIt tol dt lambda];
  tic
  T = difFinitas(xnode, model, cb, et);
  t_explicito = toc

  plot(xnode, T, 'linewidth', 1.5, 'r', 'DisplayName', 'Explicito')

  printf("Temperatura en nodo final x=L2 (Explicito): T = %.5f\n", T(end));

  % Resolucion implicita:
  et = [2 maxIt tol dt lambda];
  tic
  T = difFinitas(xnode, model, cb, et);
  t_implicito = toc

  printf("Temperatura en nodo final x=L2 (Implicito): T = %.5f\n", T(end));

  plot(xnode, T, 'linewidth', 1.5, 'g', 'DisplayName', 'Implicito')

  xlabel("xi")
  ylabel("Temperatura")
  title("Comparación: Estacionaria - Explicito - Implicito")
  legend('Location', 'southeast')

