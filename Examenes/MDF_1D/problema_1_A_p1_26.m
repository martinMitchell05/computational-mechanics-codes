clc; clear all;

  L1 = 0;
  L2 = 1;
  N = 5;

  %T_1 := Neumann
  % T_2 := Robin
  q = -10;
  h = 5;
  T_inf = 0;

  k = 2;
  c = 3;
  G = 60;
  dx = (L2 - L1)/(N-1);
  cp = 1;
  p = 4;

  T_n = zeros(N,1);   % Condición inicial T(x, 0)
  x = L1:dx:L2;

  model = {k, c, p, cp, G, T_n};
  cb = [2 q -1; 3 h T_inf];
  et = [0 0 0 0 0];

  tic;
  disp("Solución estado estacionario del modelo:")
  T = difFinitas(x, model, cb, et);
  disp(T)
  tiempo_calculo = toc;


  figure(1, 'name', "Grafica de curva de temperatura")
  plot(x, T, 'k--*', 'markersize', 7, 'linewidth', 1.5, 'color', 'b')
  grid on; hold on;

  ylabel("Temperatura")
  xlabel(" xi ")
  legend("MDF", 'location', 'south');


  % Inciso A. 2. Resolver transiente con backward Euler (implicito)

  tol = 1e-7;
  dt = 0.5;
  maxIt = 2;

  et = [2, maxIt, tol, dt, 0.0];

  [T_impl, T_vec] = difFinitas(x, model, cb, et);

  plot(x, T_impl, 'k*-', 'markersize', 7, 'linewidth', 2, 'DisplayName', 'Implicita');




##
##  % Indice del nodo x = 0.5:
##  idx_05 = find(xnode == 0.5);
##
##  % Buscar en T_vec cuando en ese nodo sea mayor que 23:
##  idx_t = find(T_vec(idx_05, :) > 23, 1, 'first');
##
##  % Calcular cual es el instante t correspondiente a esa columna:
##  printf('Instante de tiempo en que T(0.5, t) > 23: t = %.2f s, Iteración %d\n', dt*(idx_t-1), idx_t-1);    % La primer iteración es la condición inicial.


