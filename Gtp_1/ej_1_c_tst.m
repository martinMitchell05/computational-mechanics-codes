% Prueba de función: Ej-1.c

clear all; clc;
  N = 16;
  L1 = 1;
  L2 = 5;

  % T_1 := Neumann = q = 2
  q = 2;
  T_2 = 0;

  k = 1;
  c = 0;
  G = @(x) 100 .* (x - 3).^2;
  dx = (L2 - L1)/(N-1);
  cp = 0;
  p = 0;

  x = L1:dx:L2;
  model = {k, c, p, cp, G}; % pasar model como un 'arreglo de celdas' para que pueda pasarle G como función
  cb = [2 q -1; 1 T_2 -1];
  et = 0;

  tic
  disp("Prueba de método 'difFinitas':")
  T = difFinitas(x, model, cb, et)
  tiempo_calculo = toc


  % Solución analítica:
  a = -1350;
  b = 1906;
  c = 2345;

  T_a = @(x) (-25.*x.^4 + 300.*x.^3 + a.*x.^2 + b.*x + c)./3;

  disp("\nSolución analítica:")
  disp(T_a(x'))

  figure(1, 'name', "Grafica comparación")
  plot(x, T, 'linewidth', 1.5, 'color', 'r')
  grid on; hold on;
  plot(x, T_a(x), 'ko-', 'markersize', 5, 'linewidth', 1.5, 'color', 'b')

  xlabel("  xi  ")
  ylabel("Temperatura")
  legend("MDF", "Analitico")
