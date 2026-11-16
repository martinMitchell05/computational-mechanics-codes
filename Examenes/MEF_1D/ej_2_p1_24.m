clc; clear all;

  % Ejercicio 2 - Parcial 1, 2024.
  L1 = 2;
  L2 = 4;
  T_1 = 10;
  q = 1;


  % Modelo:
  model.rho = 1;
  model.cp = 1;
  model.k = 2;
  model.c = 0;
  model.G = @(x) 100 .* (4 - x);
  model.tipo_e = 'lineal';
  model.T0 = [10, 39.848, 52.149, 54.173, 52.852];   % Condición inicial inciso 'b' (T^n)

  % Malla:
  N = 5;   % cantidad ptos
  dx = (L2-L1)/(N-1);
  xnode = L1:dx:L2;

  % Condiciones de borde:
  cb = [1 T_1 0;2 q 0];

  % Datos temporales:
  dt = 0.1;
  tol = 1e-5;
  maxIt = 500;
  t0 = 0;
  tF = 1;

  et = [2, dt, tol, maxIt, t0, tF];

  [T, tvec, xdof] = elemFinitos(xnode, model, cb, et);

  disp("Solución en T^(n+1) encontrada por esquema implicito:")
  disp(T(:, 2));   % la primer columna sería T^n


  % ---- Inciso c). Encontrar temperatura en x = 3.2 según resultado anterior:

  % Interpolación del punto tomando x = 3 y x = 3.5:
  % Calcular los pesos de cada función de prueba N del elemento involucrado (e3): N3 y N3.5:
  w_N3 = (3.5 - 3.2)/dx;
  w_N35 = (3.2 - 3)/dx;
  Tx_32 = w_N3*T(3,2) + w_N35*T(4,2);
  printf('Temperatura en x = 3.2: T = %g\n', Tx_32);

  % ---- Inciso d). Calcular el flujo (valor y dirección) en x = 3.2:

  % flujo := q = -k * dT/dx
  dT = (T(4,2) - T(3,2))/dx;
  flujo_x32 = -model.k*dT;
  printf('Flujo en x = 3.2: q = %g\n', flujo_x32);
  % Flujo negativo ---> dirección '-x', el calor va de la zona más caliente (derecha de la barra) a la más fría (izquierda de la barra)
  % Se puede ver que la temperatura 'sube' de izquierda a derecha en este caso, por eso para la izquierda sería enfriamiento.

