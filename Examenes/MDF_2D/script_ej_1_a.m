%% --- Script general de inicialización para MDF-2D ---

% El script arma la estructura de la grilla del problema + condiciones en los bordes, luego
% llama a la función propuesta en clase: fdm2d_solve y grafica los resultados (fdm2d_graph_mesh)

% --------------------------------------------------------------------------------------------------
  close all; clear all; more off; clc;

  % --- Parámetros / Variables del problema ---
  Lx = 1.0; % Longitud 'x'
  Ly = 2.0; % Longitud 'y'
  Nx = 16;   % Nodos en el eje x
  Ny = 16;   % Nodos en el eje y

  k = 2;
  c = 4;
  G = 300;  % Cambiar si G es G(x)

  ts = 2;
  rho = 1;
  p = 1;
  cp = 2;
  maxit = 500;
  tol = 1e-6;


  % xnode:
  x = linspace(0, Lx, Nx);
  y = linspace(0, Ly, Ny);
  % orden en meshgrid para que 'y' varíe más rápido (avanza primero en dirección 'y')
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


  % Condiciones de borde (Indexación lógica)
  tol_borde = 1e-6;

  % Encuentra todos los nodos que están en los bordes del rectángulo
  idx_borde = find(xnode(:,1) <= tol_borde | xnode(:,1) >= Lx-tol_borde | ...
                   xnode(:,2) <= tol_borde | xnode(:,2) >= Ly-tol_borde);

                   % Corregir para los bordes T=0
  DIR = [idx_borde, 20 * ones(length(idx_borde), 1)];   % Formato [ indice, valor_dirichlet ]
  NEU = [];
  ROB = [];


  disp('---------------------------------------------------------------');
  disp('Inicializando modelo de datos...');

  model.k = k * ones(model.nnodes, 1);
  model.c = c * ones(model.nnodes, 1);
  model.G = G * ones(model.nnodes, 1);

  % Parámetros temporales
  % ts: [0] Explícito, [1] Implícito, [X] Estacionario
  model.ts = ts;
  model.rho = p;
  model.cp = cp;
  model.maxit = maxit;
  model.tol = tol;

  % Condición inicial
  model.PHI_n = mean(DIR(:,2))*ones(model.nnodes,1);

  disp('Iniciando el método numérico...');

  % Llamada principal al Método de Diferencias Finitas
  [PHI,Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

  disp('Finalizada la ejecución del método numérico.');

  disp('---------------------------------------------------------------');
  disp('Iniciando el post-procesamiento...');

  % mode ---> modo de visualización:
  %           [0] 2D - Con malla
  %           [1] 3D - Con malla
  %           [2] 2D - Sin malla
  %           [3] 3D - Sin malla
  % graph --> tipo de gráfica:
  %           [0] Temperatura (escalar)
  %           [1] Flujo de Calor (vectorial)
  %           [2] Flujo de Calor eje-x (escalar)
  %           [3] Flujo de Calor eje-y (escalar)
  %           [4] Magnitud de Flujo de Calor (escalar)
  mode = 1;
  graph = 0;
  fdm2d_graph_mesh(full(PHI),Q,xnode,icone,mode,graph);

  disp('Finalizado el post-procesamiento.');


