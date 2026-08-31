% Prueba de función: Ej-1.b

clear all;
  N = 16;
  L1 = 0;
  L2 = 2;

  T_1 = 100;
  % T_2 := Neumann = q = 0
  q = 0;

  k = 1;
  c = 1;
  G = 0;
  dx = (L2 - L1)/(N-1);
  cp = 0;
  p = 0;

  x = L1:dx:L2;
  model = {k, c, p, cp, G};
  cb = [1 T_1 -1; 2 q -1];
  et = 0;

  tic
  disp("Prueba de método 'difFinitas':")
  T = difFinitas(x, model, cb, et)
  tiempo_calculo = toc


  % Analitico:
  T_a = @(x) (100 .* exp(-x).*(exp(2 .*x) + exp(4)))./(1 + exp(4));

  disp("\nSolución Analítica:")
  sol_a = T_a(x')
##  sol_a(end-2:end)

  figure(1, 'name', "Grafico comparación")
  plot(x, T, 'linewidth', 1.5, 'r-')
  grid on; hold on;
  plot(x, sol_a, 'ko-', 'markersize', 5, 'markerfacecolor', 'b', 'linewidth', 1.5, 'color', 'c')

  ylabel("Temperatura")
  xlabel(" xi ")
  legend("MDF", "Analitico")
