  clear all; clc;
  % Ejercicio 2 - GTP MDF-2D.

  % Inciso d). Resolver estado estacionario pero con fuente G(x,y) cambiante en el tiempo:
  % --- Dividir el problema en 2 etapas: Primero con la fuente activa (en este caso), Etapa final con fuente desactivada (en este caso).

  % --- Parámetros / Variables del problema ---
  Lx = 3.0; % Longitud total del rectángulo envolvente en 'x'
  Ly = 6.0; % Longitud total del rectángulo envolvente en 'y'

  % Número de nodos en la caja contenedora
  % Elegimos resoluciones que coincidan exactamente con las aristas (x=1, x=3, y=2, y=6)
  Nx = 31;   % dx = 3.0 / (31-1) = 0.1
  Ny = 61;   % dy = 6.0 / (61-1) = 0.1

  k = 1;
  c = 0;

  ts = 1;
  rho = 1;
  cp = 1;
  tol = 1e-5;
  dt = 0.05;

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

  % --- Asignación espacial de la fuente G(x, y) ---
  tol_borde = 1e-5;
  model.G = zeros(model.nnodes, 1); % Inicializamos todo el dominio con G = 0

  % Máscara lógica para seleccionar los nodos en 0 <= x <= 1 y 0 <= y <= 3
  mask_G = (xnode(:,1) >= 0.0 - tol_borde & xnode(:,1) <= 1.0 + tol_borde) & ...
           (xnode(:,2) >= 0.0 - tol_borde & xnode(:,2) <= 3.0 + tol_borde);

  % Asignar G = 100 únicamente a los nodos pertenecientes a ese subdominio
  model.G(mask_G) = 100;

  % Crear mapa de reindexación de nodos (de la grilla global a la L reducida)
  mapa_nodos = zeros(max(nodos_usados), 1);
  mapa_nodos(nodos_usados) = (1:model.nnodes)';
  icone = mapa_nodos(icone_raw);

  % --- 4. Asignación de Condiciones de Borde sobre xnode activo ---

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

  % Parámetros temporales
  model.ts = ts;
  model.rho = rho;
  model.cp = cp;
  model.maxit = round(2.0/dt);    % Primer etapa: calcular hasta que se apaga la fuente = Tf/dt
  model.tol = tol;
  model.dt = dt;

  % Condición inicial
  model.PHI_n = mean(DIR(:,2))*ones(model.nnodes,1);

  disp('Iniciando el método numérico - ETAPA 1...');

  % Llamada principal al Método de Diferencias Finitas
  [PHI_1,Q_1] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

  disp('Finalizada la ejecución del método numérico - ETAPA 1.');

  % =========================================================================
  % ETAPA 2: Transitorio sin Fuente (G = 0) de t = 2.0 s hasta Estacionario
  % =========================================================================
  disp('--- Iniciando Etapa 2 (t > 2 s) con Fuente Apagada ---');

  % 1. La condición inicial de la Etapa 2 es la solución en t = 2.0 s (última columna)
  model.PHI_n = PHI_1(:, end);

  % 2. Se apaga la fuente G en todo el dominio
  model.G = zeros(model.nnodes, 1);

  % 3. Criterio de corte por tolerancia de estacionario (1e-5)
  model.maxit = 500;
  model.tol = 1e-5;
  model.dt = dt;

  [PHI_2, Q_2] = fdm2d(xnode, icone, DIR, NEU, ROB, model);


  % =========================================================================
  % Unificación de Resultados y Post-procesamiento
  % =========================================================================
  % Se concatenan las soluciones omitiendo la columna repetida de t = 2.0 s
  PHI_total = [PHI_1, PHI_2(:, 2:end)];
  Q_total   = [Q_1,   Q_2(:, 3:end)]; % 2 columnas por paso en Q

  fprintf('Proceso completado. Total de pasos temporales simulados: %d\n', size(PHI_total, 2) - 1);
  printf("Tiempo total transcurrido hasta alcanzar el estado estacionario: t = %.2f\n", dt*(size(PHI_total,2)-1));

  disp('---------------------------------------------------------------');
  disp('Iniciando el post-procesamiento...');

  mode = 0; % Gráfica 2D con malla
  graph = 0; % Superficie de Temperatura
  fdm2d_graph_mesh(full(PHI_total),Q_total,xnode,icone,mode,graph);

  disp('Finalizado el post-procesamiento.');


  % ---- Informar los valores de temperatura en un rango de nodos:
  % y = 6 (x en [0, 1]) a t = 1, 5, 10 s
  % =========================================================================

  tol_borde = 1e-5;
  idx_linea = find(abs(xnode(:,2) - 6.0) <= tol_borde & ...
                   xnode(:,1) >= 0.0 - tol_borde & ...
                   xnode(:,1) <= 1.0 + tol_borde);

  % 2. Ordenar los nodos por coordenada x creciente
  [x_coords, order] = sort(xnode(idx_linea, 1));
  idx_linea = idx_linea(order);

  % 3. Identificar índices de columna para t = 1 s, t = 5 s, t = 10 s
  col_t1  = round(1.0  / dt) + 1;
  col_t5  = round(5.0  / dt) + 1;
  col_t10 = round(10.0 / dt) + 1;

  % 4. Extraer temperaturas en la línea elegida
  T_t1  = PHI_total(idx_linea, col_t1);
  T_t5  = PHI_total(idx_linea, col_t5);
  T_t10 = PHI_total(idx_linea, col_t10);

  % 5. Mostrar la tabla de resultados en pantalla
  fprintf('\n=======================================================\n');
  fprintf('  Perfil de Temperatura T(x, y=6) para t = 1s, 5s, 10s\n');
  fprintf('=======================================================\n');
  fprintf('  x (m)  |  T(t=1s) [ºC]  |  T(t=5s) [ºC]  |  T(t=10s) [ºC]\n');
  fprintf('-------------------------------------------------------\n');
  for k = 1:length(idx_linea)
      fprintf('  %5.2f  |   %10.4f   |   %10.4f   |   %10.4f\n', ...
              x_coords(k), T_t1(k), T_t5(k), T_t10(k));
  end
  fprintf('=======================================================\n');

  % 6. Graficar las curvas de temperatura vs x
  figure;
  plot(x_coords, T_t1,  '-or', 'LineWidth', 1.5, 'DisplayName', 't = 1 s'); hold on;
  plot(x_coords, T_t5,  '-sb', 'LineWidth', 1.5, 'DisplayName', 't = 5 s');
  plot(x_coords, T_t10, '-^g', 'LineWidth', 1.5, 'DisplayName', 't = 10 s');
  grid on;
  xlabel('Coordenada x [m]');
  ylabel('Temperatura [ºC]');
  title('Evolución de la Temperatura en el borde y = 6 m');
  legend('Location', 'best');
