clc; clear all;

  L1 = 0;
  L2 = 1;
  %T_1 := Neumann
  q = -10;

  %T_2 := Robin
  h = 5;
  T_inf = 0;

  % Malla:
  N = 11;   % cantidad ptos
  dx = (L2-L1)/(N-1);
  xnode = L1:dx:L2;

  % Modelo:
  model.rho = 4;
  model.cp = 1;
  model.k = 2;
  model.c = 3;
  model.G = 60;
  model.tipo_e = 'lineal';
  model.T0 = zeros(N,1);

  % Condiciones de borde:
  cb = [2 q 0;3 h T_inf];

  et = [0, 0, 0, 0, 0, 0];

  [T, tvec, xdof] = elemFinitos(xnode, model, cb, et);

  disp("Solución del sistema estado estacionario:");
  disp(T);

  % Mostrar solo los 5 puntos de inciso A.

  T_025 = (T(4) + T(3))/2;
  T_075 = (T(8) + T(9))/2;

  T_5_ptos = [T(1); T_025; T(6); T_075; T(end)];
  disp("\nSolución solo en los 5 puntos de inciso A, N = 11 nodos: ");
  disp(T_5_ptos);

  figure(1, 'name', "Grafica de curva de temperatura")
  plot(xnode, T, 'k--*', 'markersize', 7, 'linewidth', 1.5, 'color', 'b')
  grid on; hold on;

  ylabel("Temperatura")
  xlabel(" xi ")
  legend("MDF", 'location', 'south');


  % Inciso D.1.2:
  N = 21;
  dx2 = (L2 - L1)/(N-1);
  x2 = L1:dx2:L2;

  disp("\nSolución estado estacionario con N = 21 nodos: ");
  T_21 = elemFinitos(x2, model, cb, et)

  T_21_5_ptos = [T_21(1); T_21(6); T_21(11); T_21(16); T_21(end)];

  disp("\nSolución en los 5 nodos, N = 21 puntos: ");
  disp(T_21_5_ptos);
