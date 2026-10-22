  clear all; clc;
  % Ejercicio 2 - GTP MDF-2D.

  % Inciso b). Resolver estado transiente:

  % --- Parámetros / Variables del problema ---
  Lx = 3.0; % Longitud total del rectángulo envolvente en 'x'
  Ly = 6.0; % Longitud total del rectángulo envolvente en 'y'

  % Número de nodos en la caja contenedora
  % Elegimos resoluciones que coincidan exactamente con las aristas (x=1, x=3, y=2, y=6)
  Nx = 31;   % dx = 3.0 / (31-1) = 0.1
  Ny = 61;   % dy = 6.0 / (61-1) = 0.1

  k = 1;
  c = 0;
  G = 100;    % Fuente de calor interna (ajustar si el enunciado especifica otro valor)

  ts = 1;   % Estado transiente: Implicito
  rho = 1;
  cp = 1;
  dt = 0.05;
  maxit = 500;  % Explicito con 500 no llega a estacionario.
  tol = 1e-5;

  % --- 1. Generación del Grid Rectangular Envolvente ---
  x = linspace(0, Lx, Nx);
  y = linspace(0, Ly, Ny);
  [X, Y] = meshgrid(x, y);
  xnode_full = [X(:), Y(:)];

  % --- 2. Generación de Icone (Filtrando la geometría en 'L') ---
  num_elem_max = (Nx-1) * (Ny-1);
  icone_raw = [];

  for i = 1:(Nx-1)
      for j = 1:(Ny-1)
          % Centroide del elemento cuadrilátero
          xc = (x(i) + x(i+1)) / 2;
          yc = (y(j) + y(j+1)) / 2;

          % Pertenece a la L si NO está en el recorte superior derecho (x > 1 e y > 2)
          if (xc <= 1.0 || yc <= 2.0)
              n1 = (i-1)*Ny + j;     % inf-izq
              n2 = i*Ny + j;         % inf-der
              n3 = i*Ny + j + 1;     % sup-der
              n4 = (i-1)*Ny + j + 1; % sup-izq

              icone_raw = [icone_raw; n1, n2, n3, n4];
          end
      end
  end

  % --- 3. Filtrado de Nodos Activos y Reindexación ---
  nodos_usados = unique(icone_raw(:));
  xnode = xnode_full(nodos_usados, :);
  model.nnodes = size(xnode, 1);

  % Crear mapa de reindexación de nodos (de la grilla global a la L reducida)
  mapa_nodos = zeros(max(nodos_usados), 1);
  mapa_nodos(nodos_usados) = (1:model.nnodes)';
  icone = mapa_nodos(icone_raw);

  % --- 4. Asignación de Condiciones de Borde sobre xnode activo ---
  tol_borde = 1e-5;

  % A) Dirichlet: Borde superior (y = 6, x en [0, 1]) -> T = 0
  idx_top = find(abs(xnode(:,2) - 6.0) <= tol_borde & xnode(:,1) <= 1.0 + tol_borde);
  DIR_top = [idx_top, zeros(length(idx_top), 1)];

  % B) Dirichlet: Borde izquierdo (x = 0, y en [0, 6)) -> T = 20
  idx_left = find(abs(xnode(:,1) - 0.0) <= tol_borde & xnode(:,2) < 6.0 - tol_borde);
  DIR_left = [idx_left, 20 * ones(length(idx_left), 1)];

  DIR = [DIR_top; DIR_left];

  % C) Robin: Borde inferior (y = 0, x en [0, 3]) -> h = 2, T_inf = 50, Dirección = 1 (Sur)
  idx_bottom = find(abs(xnode(:,2) - 0.0) <= tol_borde & xnode(:,1) <= 3.0 + tol_borde);
  ROB = [idx_bottom, 2 * ones(length(idx_bottom), 1), 50 * ones(length(idx_bottom), 1), 1 * ones(length(idx_bottom), 1)];

  % D) Neumann: Flujos nulos q = 0
  % Pared interna vertical (x = 1, y en [2, 6]) -> Dir = 2 (Este)
  idx_in_v = find(abs(xnode(:,1) - 1.0) <= tol_borde & xnode(:,2) >= 2.0 - tol_borde & xnode(:,2) <= 6.0 + tol_borde);
  NEU_in_v = [idx_in_v, zeros(length(idx_in_v), 1), 2 * ones(length(idx_in_v), 1)];

  % Pared interna horizontal (y = 2, x en [1, 3]) -> Dir = 3 (Norte)
  idx_in_h = find(abs(xnode(:,2) - 2.0) <= tol_borde & xnode(:,1) >= 1.0 - tol_borde & xnode(:,1) <= 3.0 + tol_borde);
  NEU_in_h = [idx_in_h, zeros(length(idx_in_h), 1), 3 * ones(length(idx_in_h), 1)];

  % Pared externa derecha (x = 3, y en [0, 2]) -> Dir = 2 (Este)
  idx_out_r = find(abs(xnode(:,1) - 3.0) <= tol_borde & xnode(:,2) >= 0.0 - tol_borde & xnode(:,2) <= 2.0 + tol_borde);
  NEU_out_r = [idx_out_r, zeros(length(idx_out_r), 1), 2 * ones(length(idx_out_r), 1)];

  NEU = [NEU_in_v; NEU_in_h; NEU_out_r];

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
  model.dt = dt;

  % Condición inicial
  model.PHI_n = mean(DIR(:,2))*ones(model.nnodes,1);

  disp('Iniciando el método numérico...');

  % Llamada principal al Método de Diferencias Finitas
  [PHI,Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

  disp('Finalizada la ejecución del método numérico.');

  disp('---------------------------------------------------------------');
  disp('Iniciando el post-procesamiento...');

  mode = 0; % Gráfica 3D con malla
  graph = 0; % Superficie de Temperatura
  fdm2d_graph_mesh(full(PHI),Q,xnode,icone,mode,graph);

  disp('Finalizado el post-procesamiento.');
