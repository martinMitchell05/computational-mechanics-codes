clc; clear all;

  % Inciso D. Problema 1.

  L1 = 0;
  L2 = 1;
  N = 11;     % D.1.1

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


  % Mostrar solo los 5 puntos de inciso A.

  T_025 = (T(4) + T(3))/2;
  T_075 = (T(8) + T(9))/2;

  T_5_ptos = [T(1); T_025; T(6); T_075; T(end)];
  disp("\nSolución solo en los 5 puntos de inciso A, N = 11 nodos: ");
  disp(T_5_ptos);

  figure(1, 'name', "Grafica de curva de temperatura")
  plot(x, T, 'k--*', 'markersize', 7, 'linewidth', 1.5, 'color', 'b')
  grid on; hold on;

  ylabel("Temperatura")
  xlabel(" xi ")
  legend("MDF", 'location', 'south');


  % Inciso D.1.2:
  N = 21;
  dx2 = (L2 - L1)/(N-1);
  x2 = L1:dx2:L2;

  disp("\nSolución estado estacionario con N = 21 nodos: ");
  T_21 = difFinitas(x2, model, cb, et)

  T_21_5_ptos = [T_21(1); T_21(6); T_21(11); T_21(16); T_21(end)];

  disp("\nSolución en los 5 nodos, N = 21 puntos: ");
  disp(T_21_5_ptos);
