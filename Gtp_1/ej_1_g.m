% Ejercicio 1.g
clc; clear all;

  L1 = 0;
  L2 = 1;

  T_1 = 50;
  % T_2 = Neumann
  q = 5;

  N = 32;
  dx = (L2 - L1)/(N-1);
  xnode = L1:dx:L2;

  k = 2;
  G = 0;
  c = -2;
  pcp = 1;
  p = 1;

  lambda = 0.97; % un lambda mas bajo puede empeorar el resultado, no mejorar
  dt = 0.4; % un dt no excesivamente bajo y da un muy buen resultado

  tol = 1e-7;
  maxIt = 7000;
  % Mientras mas iteraciones (y/o) mas baja tolerancia : mejor va a aproximar el Explicito
  % Pero eso implica que tarde mucho más tiempo para los mismos puntos

  model = {k, c, p, pcp, G};
  cb = [1 T_1 0; 2 q 0];

  disp("Solución por Diferencias Finitas:\n")

  % Solucion estacionaria
  et = [0 0 0 0 0];

  tic
  T_est = difFinitas(xnode, model, cb, et);
  t_estacionaria = toc

  figure(1)
  plot(xnode, T_est, 'ko-', 'markersize', 7, 'linewidth', 1.5, 'color', 'b', 'DisplayName', 'Estacionario')
  grid on; hold on;

  % Solucion explicita
  et = [1 maxIt tol dt lambda];

  tic
  T_expl = difFinitas(xnode, model, cb, et);
  t_explicito = toc

  plot(xnode, T_expl, 'linewidth', 1.5, 'r', 'DisplayName', 'Explicito')

  % Solucion implicita
  et = [2 maxIt tol dt lambda];

  tic
  T_impl = difFinitas(xnode, model, cb, et);
  t_impl = toc

  plot(xnode, T_impl, 'linewidth', 1.5, 'm', 'DisplayName', 'Implicito')


  % Calculo de errores:
  error_est_expl = norm(T_est - T_expl, 2);
  error_est_impl = norm(T_est - T_impl, 2);
  error_expl_impl = norm(T_expl - T_impl, 2);

  printf("Errores entre métodos: (norma 2)\nEstacionario - Explicito = {%d}\nEstacionario - Implicito = {%d}\nExplicito - Implicito = {%d}\n\n", error_est_expl, error_est_impl, error_expl_impl)


  % Solucion analitica:
  T_a = @(x) 73.2433.*sin(x) + 50.*cos(x);
  sol_a = T_a(xnode);

  plot(xnode, sol_a, 'linewidth', 2, 'g--', 'DisplayName', 'Analitica')
  xlabel("xi")
  ylabel("Temperatura")
  title("Grafica comparacion: Estacionario - Explicito - Implicito - Analitica")
  legend('Location', 'southeast')

  % Calculo de errores:
  error_est_an = norm(sol_a - T_est, 2);
  error_expl_an = norm(sol_a - T_expl, 2);
  error_impl_an = norm(sol_a - T_impl, 2);
  printf("Errores entre métodos: (norma 2)\nAnalitica - Estacionario = {%d}\nAnalitica - Explicito = {%d}\nAnalitica - Implicito = {%d}\n\n", error_est_an, error_expl_an, error_impl_an)

