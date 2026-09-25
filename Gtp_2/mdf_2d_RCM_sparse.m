  % =========================================================================
  % RESOLUCIÓN DE LA ECUACIÓN DE REACCIÓN-DIFUSIÓN 2D POR MDF (SPARSE + RCM)
  % Modelo: - (d2phi/dx2 + d2phi/dy2) + c*phi = Q/k
  % =========================================================================
  clc; clear; tic;

  % 1. Parámetros físicos y geométricos
  Lx = 4.0; Ly = 1.0;      % Dimensiones del dominio
  k = 2.0;                 % Conductividad / Difusividad
  c = 10.0;                % Coeficiente reactivo
  Q_val = 100.0;           % Fuente volumétrica
  phi_boundary = 0.0;     % Condición Dirichlet en frontera

  % 2. Discretización espacial
  L = 30; M = 30;          % Número de intervalos en x e y
  dx = Lx / L; dy = Ly / M;
  h = dx;                  % Asumimos dx = dy

  % Número de nodos interiores
  nx = L - 1; ny = M - 1;
  N = nx * ny;             % Total de grados de libertad

  % 3. Construcción eficiente de la matriz dispersa K (Triplets)
  row = []; col = []; val = [];
  f = zeros(N, 1);

  for i = 1:nx
      for j = 1:ny
          % Mapeo de coordenadas (i,j) a índice global k (numeración por y)
          k_idx = (i - 1) * ny + j;

          % Diagonal principal (Stencil de 5 puntos + Término Reactivo)
          row(end+1) = k_idx; col(end+1) = k_idx;
          val(end+1) = 4.0 + (c / k) * h^2;

          % Vector de carga interior
          f(k_idx) = (h^2 * Q_val) / k;

          % Vecino izquierdo (i-1, j)
          if i > 1
              row(end+1) = k_idx; col(end+1) = (i - 2) * ny + j;
              val(end+1) = -1.0;
          else
              f(k_idx) = f(k_idx) + phi_boundary; % Dirichlet contorno izquierdo
          end

          % Vecino derecho (i+1, j)
          if i < nx
              row(end+1) = k_idx; col(end+1) = i * ny + j;
              val(end+1) = -1.0;
          else
              f(k_idx) = f(k_idx) + phi_boundary; % Dirichlet contorno derecho
          end

          % Vecino inferior (i, j-1)
          if j > 1
              row(end+1) = k_idx; col(end+1) = (i - 1) * ny + (j - 1);
              val(end+1) = -1.0;
          else
              f(k_idx) = f(k_idx) + phi_boundary; % Dirichlet contorno inferior
          end

          % Vecino superior (i, j+1)
          if j < ny
              row(end+1) = k_idx; col(end+1) = (i - 1) * ny + (j + 1);
              val(end+1) = -1.0;
          else
              f(k_idx) = f(k_idx) + phi_boundary; % Dirichlet contorno superior
          end
      end
  end

  % Generación de la matriz dispersa
  K_sparse = sparse(row, col, val, N, N);

  % 4. Optimización del Ancho de Banda mediante Reverse Cuthill-McKee (RCM)
  p = symrcm(K_sparse);          % Vector de permutación RCM
  K_rcm = K_sparse(p, p);        % Matriz reordenada
  f_rcm = f(p);                  % Vector reordenado

  % 5. Resolución del sistema lineal de alta eficiencia
  phi_rcm = K_rcm \ f_rcm;

  % Deshacer la permutación RCM para mapear al dominio físico original
  phi_vec = zeros(N, 1);
  phi_vec(p) = phi_rcm;

  % 6. Reconstrucción de la matriz 2D para visualización
  PHI = zeros(ny, nx);
  for i = 1:nx
      for j = 1:ny
          k_idx = (i - 1) * ny + j;
          PHI(j, i) = phi_vec(k_idx);
      end
  end

  % 7. Gráfico de Contorno / Superficie
  figure;
  [X, Y] = meshgrid(linspace(dx, Lx-dx, nx), linspace(dy, Ly-dy, ny));
  surf(X, Y, PHI);
  title('Distribución 2D de la Solución (MDF Centrado 5 Puntos)');
  xlabel('Posición X'); ylabel('Posición Y'); zlabel('\phi(x,y)');
  colorbar; grid on;

  tiempo_script = toc
