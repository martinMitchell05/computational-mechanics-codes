  clear all; clc;
  % Ejercicio 2 - P1 2025 - Placa 1x1

  % --- Parámetros / Variables del problema ---
  Lx = 1; % Longitud total del rectángulo envolvente en 'x'
  Ly = 1; % Longitud total del rectángulo envolvente en 'y'

  % Número de nodos en la caja contenedora
  Nx = 7;   % dx = (L2x - L1x) / (Nx-1)
  Ny = 7;   % dy = (L2y - L1y) / (Ny-1)

  k = 10;
  c = 0;
  G = 500;

  ts = 0;   % Esquema explicito
  rho = 1;
  cp = 1;
  maxit = 5000;
  tol = 1e-5;

  % --- 1. Generación del Grid Rectangular Envolvente ---
  x = linspace(0, Lx, Nx);
  y = linspace(0, Ly, Ny);
  [X, Y] = meshgrid(x, y);
  xnode = [X(:), Y(:)];
  model.nnodes = size(xnode, 1);

  % Icone:
  num_elem = (Nx-1) * (Ny-1);
  icone = zeros(num_elem, 4);
  idx = 1;
  for i = 1:(Nx-1)
      for j = 1:(Ny-1)
          n1 = (i-1)*Ny + j; % Nodo inferior izquierdo del elemento

          % Recorrido antihorario: inf-izq, inf-der, sup-der, sup-izq
          icone(idx, :) = [n1, n1+Ny, n1+Ny+1, n1+1];
          idx = idx + 1;
      end
  end


  % --- 4. Asignación de Condiciones de Borde sobre xnode activo ---
  tol_borde = 1e-5;

  % A) Dirichlet: Borde superior (T = 150)
  idx_top = find(abs(xnode(:,2) - Ly) <= tol_borde & ...
                    xnode(:,1) <= Lx + tol_borde & xnode(:, 1) >= 0.0 - tol_borde);
  DIR_top = [idx_top, 150*ones(length(idx_top), 1)];

  % B) Dirichlet: Borde izquierdo (T = 100)
  idx_left = find(abs(xnode(:,1) - 0.0) <= tol_borde & ...
                    xnode(:,2) >= 0.0 - tol_borde & xnode(:,2) < Ly - tol_borde);
  DIR_left = [idx_left, 100 * ones(length(idx_left), 1)];

  DIR = [DIR_top; DIR_left];

  % C) Robin: Borde inferior (h = 50, T_inf = 20, Dirección = 1 (Sur)
  idx_bottom = find(abs(xnode(:,2) - 0.0) <= tol_borde & ...
                       xnode(:,1) >= 0.0 - tol_borde & xnode(:,1) <= Lx + tol_borde);
  ROB = [idx_bottom, 50 * ones(length(idx_bottom), 1), 20 * ones(length(idx_bottom), 1), 1 * ones(length(idx_bottom), 1)];

  % D) Neumann: Flujos nulos q = 0, Borde derecho, Dir = 2 (Este)
  idx_der = find(abs(xnode(:,1) - Lx) <= tol_borde & ...
                    xnode(:,2) >= 0.0 - tol_borde & xnode(:,2) <= Ly + tol_borde);
  NEU_der = [idx_der, zeros(length(idx_der), 1), 2 * ones(length(idx_der), 1)];

  NEU = [NEU_der];

  disp('---------------------------------------------------------------');
  disp('Inicializando modelo de datos...');

  model.k = k * ones(model.nnodes, 1);
  model.c = c * ones(model.nnodes, 1);
  model.G = G * ones(model.nnodes, 1);

  % Parámetros temporales
  model.ts = ts;
  model.rho = rho;
  model.cp = cp;
  model.maxit = maxit;
  model.tol = tol;

  % Condición inicial
  model.PHI_n = 125*ones(model.nnodes,1);

  disp('Iniciando el método numérico...');

  % Llamada principal al Método de Diferencias Finitas
  [PHI,Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

  disp('Finalizada la ejecución del método numérico.');

##  disp('---------------------------------------------------------------');
##  disp('Iniciando el post-procesamiento...');
##
##  mode = 0; % Gráfica 2D con malla
##  graph = 0; % Superficie de Temperatura
##  fdm2d_graph_mesh(full(PHI),Q,xnode,icone,mode,graph);
##
##  disp('Finalizado el post-procesamiento.');

  % Tiempo al estado estacionario:
  alpha = min(model.k)/(model.rho*model.cp);
  nd = 2;                     % porque estamos en 2D
  delta = Lx/(Nx-1);
  dt = 0.25*(delta)^2/(alpha*nd);
  t_est = dt * (size(PHI,2)-1);
  printf('\nTiempo hasta el estado estacionario: t = %.5f s\n', t_est);


