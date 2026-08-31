% Prueba de función: Ej-1.d

clear all; clc;
  N = 24;
  L1 = 0;
  L2 = 1;

  T_1 = 10;
  % T_2 := Robin
  h = 0.2;
  T_inf = 50;

  k = 1;
  c = 1;
  G = 50;
  dx = (L2 - L1)/(N-1);
  cp = 0;
  p = 0;

  x = L1:dx:L2;
  model = {k, c, p, cp, G};
  cb = [1 T_1 -1; 3 h T_inf];
  et = 0;

  tic
  disp("Prueba de método 'difFinitas':")
  T = difFinitas(x, model, cb, et)
  tiempo_calculo = toc


  % Solución analítica:
  T_a = @(x) -36.6897*exp(-x) - 3.3103*exp(x) + 50;

  disp("\nSolución analítica:")
  disp(T_a(x'))

  figure(1, 'name', "Grafico comparación")
  plot(x, T, 'linewidth', 1.5, 'r-')
  grid on; hold on;
  plot(x, T_a(x), 'ko-', 'markersize', 5, 'markerfacecolor', 'b', 'linewidth', 1.5, 'color', 'c')

  ylabel("Temperatura")
  xlabel(" xi ")
  legend("MDF", "Analitico")
