function T = difFinitas(xnode, model, cb, et)

  %% xnode : vector de coordenadas
  %% model : vector con constantes del modelo (k, c, p, cp, G)
  %% cb : matriz con condiciones de borde 2x3:
    % La primera fila indica la condición de borde del lado izquierdo y la segunda fila, la del derecho.
    % La primera columna indica el tipo de condición de borde: 1-Dirichlet, 2-Neumann, 3-Robin.
    % La segunda columna contiene el valor de temperatura, flujo, o coeficiente de convección h dependiendo el dato de la primera columna.
    % La tercera columna será de valor -1 para la condición 1 y 2, y tendrá el valor de temperatura externa en caso de la condición 3.

  %% et : esquema temporal a usar o estado estacionario [tipo, maxIt, tol, dt]
    % tipo: 0: estacionario - 1: explicito - 2: implicito, si tipo:2 ---> usa dt
  %% ---> Falta implementar.

  % Establecer constantes:
  N = length(xnode);
  dx = xnode(2) - xnode(1);

  k = model{1};
  c = model{2};
  p = model{3};
  cp = model{4};
  G = model{5};

  % --- Ver si la fuente es constante o función del nodo ---
  if isa(G, 'function_handle') % G era función
    G_vec = G(xnode)(:);
  elseif length(G) == N % G era un vector de valores
    G_vec = G(:);
  else
    G_vec = G * ones(N, 1);  % G era constate
  end

  %% ---- Ver valor de et ----
  %% --- switch para valor 'et' ---


  % ---- Mismas ecuaciones sin importar los bordes ----
  % Armar matriz K
  K = diag(-ones(N-1,1),-1) + diag((2 + c*dx^2/k)*ones(N,1)) + diag(-ones(N-1,1),1);

  % Armar vector f
  f = zeros(N,1);

  % Reemplazar en nodos interiores:
  f(2:N-1) = G_vec(2:N-1)/k*dx^2;

  % ---- Reemplazar ecuaciones de borde ----
  % Borde izquierdo:
  borde_izq = cb(1, :);
  switch (borde_izq(1))
    case 1 % Dirichlet
      f(1) = borde_izq(2);
      K(1,1:2) = [1 0];

    case 2 % Neumann
      q = borde_izq(2);
      node_fic_K = -2;
      node_fic_f = G_vec(1)/k*dx^2 - 2*dx*q/k;
      f(1) = node_fic_f;
      K(1,2) = node_fic_K;

    case 3 % Robin
      T_inf = borde_izq(3);
      h = borde_izq(2);
      f(1) = G_vec(1)/k*dx^2 + 2 * dx * h * T_inf / k;
      K(1, 1:2) = [(2 + c*dx^2/k + 2*dx*h/k), -2];

  endswitch

  % Borde derecho:
  borde_der = cb(2, :);
  switch (borde_der(1))
    case 1 % Dirichlet
      f(end) = borde_der(2);
      K(end, end-1:end) = [0 1];

    case 2 % Neumann
      q = borde_der(2);
      node_fic_f = G_vec(end)/k*dx^2 - 2*dx*q/k;
      f(end) = node_fic_f;
      K(end, end-1) = -2;

    case 3 % Robin
      T_inf = borde_der(3);
      h = borde_der(2);
      f(end) = G_vec(end)/k*dx^2 + 2 * dx * h * T_inf / k;
      K(end, end-1:end) = [-2 (2 + c*dx^2/k + 2*dx*h/k)];

  endswitch

  % Resolver sistema de ecuaciones:
  T = K\f;

endfunction

