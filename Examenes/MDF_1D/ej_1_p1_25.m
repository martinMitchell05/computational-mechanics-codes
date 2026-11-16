clc; clear all;

  % Ejercicio 1 - P1 2025:

  % Inciso 1). Resolver el estado estacionario del modelo:

  N = 4;    % 4 nodos
  L1 = 0;
  L2 = 1;

  T_1 = 0;
  % T_2 := Robin
  h = 1;
  T_inf = 0;

  k = 1;
  c = 0;
  G = 100;
  dx = (L2 - L1)/(N-1);
  cp = 1;
  p = 1;

  T_0 = @(x) 10 .*x .*(1 - x).^2;   % Condición inicial T(x, 0)
  x = L1:dx:L2;

  model = {k, c, p, cp, G};
  cb = [1 T_1 -1; 3 h T_inf];
  et = [0 0 0 0 0];

  tic;
  disp("Solución estado estacionario del modelo:")
  T = difFinitas(x, model, cb, et);
  disp(T)
  tiempo_calculo = toc;

  T_ex = @(x) 25 .* x .*(3 - 2.*x);

  sol_ex = T_ex(x)(:);    % OJO: hacerlo columna para que tenga la misma forma que el vector 'T'

  err_rel = norm(sol_ex - T, 2)/norm(sol_ex, 2);
  printf('Error relativo del modelo estacionario por MDF a la solución exacta: E = %g\n', err_rel);

  figure(1, 'name', "Grafica de curva de temperatura")
  plot(x, T, 'k--*', 'markersize', 7, 'linewidth', 1.5, 'color', 'b')
  grid on; hold on;
  plot(x, sol_ex, 'k-.', 'linewidth', 2);

  ylabel("Temperatura")
  xlabel(" xi ")
  legend("MDF", "Exacta")


  % Inciso 2). Resolver estado transiente por esquema implicito
  % lambda = 40
  N = 41;

  dx = (L2 - L1)/(N-1);
  tol = 1e-4;
  dt = 40*(dx*dx)/k;    % sale de formula del lambda
  xnode = L1:dx:L2;
  T_n = T_0(xnode);
  maxIt = 500;

  et = [2, maxIt, tol, dt, 0.0];

  [T_impl, T_vec] = difFinitas(xnode, model, cb, et);

  figure(2)
  plot(xnode, T_impl, 'k*-', 'markersize', 7, 'linewidth', 2, 'DisplayName', 'Implicita');

  % Indice del nodo x = 0.5:
  idx_05 = find(xnode == 0.5);

  % Buscar en T_vec cuando en ese nodo sea mayor que 23:
  idx_t = find(T_vec(idx_05, :) > 23, 1, 'first');

  % Calcular cual es el instante t correspondiente a esa columna:
  printf('Instante de tiempo en que T(0.5, t) > 23: t = %.2f s, Iteración %d\n', dt*(idx_t-1), idx_t-1);    % La primer iteración es la condición inicial.


