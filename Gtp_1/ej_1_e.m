% Ejercicio 1.e
clc; clear all;

  L1 = 5;
  L2 = 10;
  % T(L1) = Robin
  h = 2;
  T_inf = 100;
  T_2 = 50;

  pcp = 1;
  k = 2;
  c = 0;
  G = @(x) x.^3;
  lambda = 0.8; % lambda para bajar del dt critico
  dt = 0.5; % dt para implicito

  N = 32;
  dx = (L2 - L1)/(N-1);

  xnode = L1:dx:L2;
  model = {k, c, 1, pcp, G};
  cb = [3 h T_inf; 1 T_2 0];
  et = [0 300 1e-7 dt lambda]; % Solucion estacionaria

  disp("Resolución temporal, Diferencias Finitas:\n")
  tic
  T = difFinitas(xnode, model, cb, et);
  tiempo_calculo_estacionaria = toc

  % graficas
  figure(1, 'name', "Gráficas de solución Diferencias Finitas")
  plot(xnode, T, 'ro-', 'markersize', 10, 'linewidth', 1.5, 'color', 'k', 'DisplayName', 'Estacionaria')
  grid on;
  hold on;

  xlabel("xi")
  ylabel("Temperatura")
##  title(sprintf("Solución Estacionaria - N = %i", N))

  % Solucion explicita : forward euler
  et = [1 3000 1e-7 dt lambda];
  tic
  T = difFinitas(xnode, model, cb, et);
  tiempo_calculo_explicito = toc

  plot(xnode, T, 'linewidth', 1.25, 'b', 'DisplayName', 'Explicito')
##  title(sprintf("Solución Temporal - Forward Euler - lambda = %f - N = %i", lambda, N))

  % Solucion implicita : backward euler
  et = [2 300 1e-7 dt lambda];
  tic
  T = difFinitas(xnode, model, cb, et);
  tiempo_calculo_implicito = toc

  plot(xnode, T, 'linewidth', 1.25, 'r', 'DisplayName', 'Implicito')


  % Solucion analitica
  T_a = @(x) -x.^5 ./ 40 + 1225.*x ./ 3 - 4600/3;
  sol_a = T_a(xnode);

  plot(xnode, sol_a, 'linewidth', 2, 'g--', 'DisplayName', 'Analitica')
  legend('Location', 'northeast')
  title("Comparación: Estacionaria - Explicito - Implicito - Analitica")


