clc; clear all;

  L1 = 0;
  L2 = 1;
  %T_1 := Neumann
  q = -10;

  %T_2 := Robin
  h = 5;
  T_inf = 0;

  % Malla:
  N = 3;   % cantidad ptos
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
