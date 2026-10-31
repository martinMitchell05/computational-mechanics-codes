clear all; %clc;

  % Ejercicio 1.c

  L1 = 1;
  L2 = 5;
  q = 2;
  T_2 = 0;

  % Modelo:
  model.rho = 0;
  model.cp = 0;
  model.k = 1;
  model.c = 0;
  model.G = @(x) 100.*(x-3).^2;
  model.tipo_e = 'lineal';
  model.T0 = 0;

  % Malla:
  N = 24;   % cantidad ptos
  dx = (L2-L1)/(N-1);
  xnode = L1:dx:L2;

  % Condiciones de borde:
  cb = [2 q 0;1 T_2 0];

  % Datos temporales:
  dt = 0;
  tol = 1e-5;
  maxIt = 500;
  t0 = 0;
  tF = 2;

  et = [0 dt, tol, maxIt, t0, tF];

  [T, tvec, xdof] = elemFinitos(xnode, model, cb, et);

  disp("Solución por MEF 1D estacionario (desde mitad hasta el final):\n");
  disp(T(N/2:end))

##  T_ex = @(x) -25.*x.^2 + 65.*x + 10;

  % --- Graficas:
  figure(1)
  plot(xnode, T, 'b*-', 'linewidth', 1.5, 'DisplayName', 'Estacionaria');
  grid on; hold on;
##  plot(xnode, T_ex(xnode), 'k--', 'linewidth', 2, 'DisplayName', 'Exacta');
  legend('location', 'southeast');
  xlabel('xi');
  ylabel(' °C ');
##  title("Comparación de solución: MEF - Exacta");

