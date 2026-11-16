clear all; clc;
  % Problema 2 - P1 2026 - Placa 2x2

  % --- Parámetros / Variables del problema ---
  Lx = 2; % Longitud total del rectángulo envolvente en 'x'
  Ly = 2; % Longitud total del rectángulo envolvente en 'y'

  % Número de nodos en la caja contenedora
  Nx = 3;   % dx = (L2x - L1x) / (Nx-1)
  Ny = 5;   % dy = (L2y - L1y) / (Ny-1)

  k = 1;
  c = 0.5;
  G = 50;

  ts = 2;   % Estado estacionario
  rho = 2;
  cp = 1;
  maxit = 5000;
  tol = 1e-3;
  dt = 0.2;


  xnode = [
    0.0, 0.0;
    1.0, 0.0;
    2.0, 0.0;
    0.0, 0.5;
    1.0, 0.5;
    2.0, 0.5;
    0.0, 1.0;
    1.0, 1.0;
    2.0, 1.0;
    0.0, 1.5;
    1.0, 1.5;
    2.0, 1.5;
    0.0, 2.0;
    1.0, 2.0;
    2.0, 2.0;
  ];

  icone = [
         1,      2,      5,      4;
         2,      3,      6,      5;
         4,      5,      8,      7;
         5,      6,      9,      8;
         7,      8,     11,      10;
         8,      9,     12,      11;
        10,     11,     14,      13;
        11,     12,     15,     14

  ];

  DIR = [
         1, 10.0000000000000000;
         2, 10.0000000000000000;
         3, 10.0000000000000000
  ];

  NEU_top = [13, 0, 3; 14, 0, 3; 15, 0, 3];
  NEU_der = [3, 8, 2; 6, 8, 2; 9, 8, 2; 12, 8, 2; 15, 8, 2];

  NEU = [ NEU_top; NEU_der];

  ROB = [1, 2, 40, 4; 4, 2, 40, 4; 7, 2, 40, 4; 10, 2, 40, 4; 13, 2, 40, 4];

  model.nnodes = size(xnode, 1);


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
  model.PHI_n = zeros(model.nnodes,1);

  disp('Iniciando el método numérico...');

  % Llamada principal al Método de Diferencias Finitas
  [PHI,Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

  disp('Finalizada la ejecución del método numérico.');

  disp('---------------------------------------------------------------');
  disp('Iniciando el post-procesamiento...');

  mode = 0; % Gráfica 2D con malla
  graph = 0; % Superficie de Temperatura
  fdm2d_graph_mesh(full(PHI),Q,xnode,icone,mode,graph);

  disp('Finalizado el post-procesamiento.');


  % Inciso B.3:
  model.ts = 1;   % Esquema implicito
  [PHI_impl, Q_impl] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

  err = norm(PHI - PHI_impl(:,end), 2);
  printf("\nNorma del error L^2 = %g\n", err);

  % Inciso A.3:

  q_x = - k * ((PHI(9)-PHI(7))/(2*1));
  q_y = - k * ((PHI(11)-PHI(5))/(2*0.5));

  printf("\nFlujos en nodo 8:\nqx = %g \nqy = %g\n\n", q_x, q_y);
