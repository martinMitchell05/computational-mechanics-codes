clear all; clc;

  % Ejercicio 1.f

  L1 = 0;
  L2 = 1;
  T_1 = 0;
  h = 2;
  T_inf = 10;

  % Modelo:
  model.rho = 2;
  model.cp = 1;
  model.k = 2;
  model.c = 2;
  model.G = @(x) 75;
  model.tipo_e = 'lineal';
  model.T0 = 0;

  % Malla:
  N = 24;   % cantidad ptos
  dx = (L2-L1)/(N-1);
  xnode = L1:dx:L2;

  % Condiciones de borde:
  cb = [1 T_1 0;3 h T_inf];

  % Datos temporales:
  dt = 0.1;
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
##  legend('location', 'southeast');
  xlabel('xi');
  ylabel(' °C ');
##  title("Comparación de solución: MEF - Exacta");

##  et = [1 dt, tol, maxIt, t0, tF];      % Metodo explicito ---> Inestable
##
##  [T, tvec, xdof] = elemFinitos(xnode, model, cb, et);
##
##  disp("Solución por MEF 1D explicito (desde mitad hasta el final):\n");
##  disp(T(N/2:end,end))
##
##  figure(2); grid on; hold on;
##  for t = 1:length(tvec)
##    plot(xdof, T(:,t), 'k*', 'linewidth', 1.2);
##    pause(0.02);
##  endfor


  et = [2 dt, tol, maxIt, t0, tF];

  [T, tvec, xdof] = elemFinitos(xnode, model, cb, et);

  disp("Solución por MEF 1D implicito (desde mitad hasta el final):");
  disp(T(N/2:end, end))


  % (a) Graficar la evolución de cada paso de tiempo
  figure(3); grid on; hold on;
  for t = 1:length(tvec)
    plot(xdof, T(:,t), 'k*-', 'linewidth', 1.2);
    pause(0.1);
  endfor


  % (b) Animación de una sola curva que se actualiza
  figure(4); clf;
  h = plot(xdof, T(:,1), 'k*-', 'linewidth', 1.2);
  axis([xdof(1) xdof(end) min(T(:)) max(T(:))]); grid on;
  xlabel('x'); ylabel('T');
  paso = max(1, floor(length(tvec)/100));      % máximo ~100 cuadros
  for j = 1:paso:length(tvec)
    set(h, 'ydata', T(:,j));
    title(sprintf('t = %.3g', tvec(j)));
    drawnow; pause(0.05);
  endfor


  % (c) Evolución temporal en un punto
  [~, i] = min(abs(xdof - 0.5));    % Evolución en x = 0.5
  figure(5); clf; plot(tvec, T(i,:), 'linewidth', 1.2); grid on;
  xlabel('t'); ylabel(sprintf('T(x = %.3g)', xdof(i)));

