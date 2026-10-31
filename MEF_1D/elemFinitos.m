function [T, tvec, xdof] = elemFinitos(xnode, model, cb, et)

  %%%
  % Implementar una función [T] = elemFinitos(xnode, model, cb, et) que resuelva el modelo completo de transferencia de calor donde:
    % - xnode es el vector de coordenadas nodales (nodos de esquina de cada elemento; los nodos medios del
    %   elemento cuadrático se generan solos en el punto medio).

    % - model es un struct que contiene todas las constantes del modelo (k, c, 𝛒, 𝒄𝒑 y G).
    %     Campos usados: model.tipo_e ('lineal' o 'cuadratico'), model.k, model.c,
    %     model.G (escalar, vector con un valor por grado de libertad, o function_handle),
    %     y solo en transiente: model.rho, model.cp y, opcionalmente, model.T0 (condición inicial: escalar,
    %     vector por grado de libertad o function_handle; si no existe se usa T0 = 0).

    % - cb es una matriz de dos filas y tres columnas que contiene los datos de las condiciones de borde.
      % La primera fila indica la condición de borde del lado izquierdo y la segunda
      % fila, la del derecho. La primera columna indica el tipo de condición de borde:
        % 1-Dirichlet, 2-Neumann, 3-Robin. La segunda columna contiene el valor de temperatura, flujo, o coeficiente de convección h
        % dependiendo el dato de la primera columna.
        % La tercera columna será de valor -1 para la condición 1 y 2, y tendrá el valor de temperatura externa en caso de la condición 3.

   % - et es un vector que indica el esquema temporal a utilizar (o resolver en estado estacionario) + datos constantes para transciente (dt, tolerancia, maxIt, t0, tF).
      % et = [esquema, dt, tol, maxIt, t0, tF]
      % esquema: 0 - Estacionario, 1 - Euler explícito (theta=0), 2 - Euler implícito (theta=1)
      % tol:   se corta la marcha en el tiempo si max|T^{n+1} - T^n| < tol (llegó al estacionario). tol = 0 lo desactiva.
      % maxIt: cantidad máxima de pasos de tiempo.

    %
    % Salidas:
    %   T    - estacionario: vector Ndof x 1.
    %          transiente:   matriz Ndof x (Npasos+1); la columna j es T en t = tvec(j). La última columna es el estado final.
    %   tvec - vector de tiempos (vacío en estacionario).
    %   xdof - coordenadas de los grados de libertad (lineal: xnode; cuadrático: esquinas + puntos medios), para graficar.
    %          Ndof = Nnod (lineal) o 2*Nnod-1 (cuadrático).
  %%%

  Nnod = length(xnode);
  Nelem = Nnod - 1;

  % -- Conectividad y coordenadas de los grados de libertad según el tipo de elemento
  if strcmp(model.tipo_e, 'lineal')
    Ndof = Nnod;
    % Conectividad: fila e = nodos globales del elemento e
    icone = [(1:Nelem)' (2:Nelem+1)'];  % e=1 == [x_1 x_2]
    xdof = xnode(:);

  elseif strcmp(model.tipo_e, 'cuadratico')
    Ndof = 2*Nelem + 1;
    % Elemento e = [esquina izq, nodo medio, esquina der] -> n = 2e-1, [n n+1 n+2]
    icone = [(1:2:2*Nelem-1)' (2:2:2*Nelem)' (3:2:2*Nelem+1)'];  % e=1 == [x_1 x_1.5 x_2]
    xn = xnode(:);
    xdof = zeros(Ndof, 1);
    xdof(1:2:end) = xn;
    xdof(2:2:end) = (xn(1:end-1) + xn(2:end))/2;

  else
    error('elemFinitos: tipo de elemento "%s" no implementado (use ''lineal'' o ''cuadratico'').', model.tipo_e);
  endif

  % -- Constantes del modelo
  k = model.k;
  c = model.c;
  G = model.G;

  esquema = et(1);
  if esquema ~= 0     % rho*cp solo hace falta en transiente
    rho = model.rho;
    cp  = model.cp;
  endif

  if isa(G, 'function_handle')
    G_vec = [];       % se evalúa dentro del lazo de elementos
  elseif length(G) == Ndof % G era un vector de valores
    G_vec = G(:);
  else
    G_vec = G * ones(Ndof, 1);  % G era constate
  endif

  A = zeros(Ndof, Ndof);  % Matriz global (K + M + H)
  C = zeros(Ndof, Ndof);  % Matriz de capacidad (solo se usa en transiente)
  b = zeros(Ndof, 1);     % Vector global (F + Q)

  % --- Ver que tipo de elemento utilizar:
  if strcmp(model.tipo_e, 'lineal')
    for e = 1:Nelem

      nod = icone(e,:);
      he = xnode(nod(2)) - xnode(nod(1));
      Ke = k/he * [1 -1; -1 1];
      Me = c*he/6 * [2 1; 1 2];
      if isa(G, 'function_handle')
        G1 = G(xnode(nod(1)));
        G2 = G(xnode(nod(2)));
        Fe = (he/6) * [2 1; 1 2] * [G1; G2];
      else
        % G interpolado con las funciones de forma (si G es constante da G*he/2*[1;1])
        Fe = (he/6) * [2 1; 1 2] * G_vec(nod);
      endif

      A(nod, nod) += (Ke + Me);
      b(nod) += Fe;

      if esquema ~= 0
        C(nod, nod) += rho*cp*he/6 * [2 1; 1 2];
      endif
    endfor

  else   % cuadratico: coordenadas naturales + integración de Gauss (como FuncCuad.m)

    % Funciones de forma para elemento master
    N = @(eta)[eta*(eta-1)*0.5, (1+eta).*(1-eta), eta*(eta+1)*0.5];
    dN_deta = @(eta)[eta-0.5, -2*eta, eta+0.5];

    % 3 puntos de gauss
    pospg = [-sqrt(15)/5 0 sqrt(15)/5];
    pespg = [5/9 8/9 5/9];

    N_eval = zeros(length(pospg), 3);
    dN_eval = zeros(length(pospg), 3);
    for pg = 1:length(pospg)
      N_eval(pg,:) = N(pospg(pg));
      dN_eval(pg,:) = dN_deta(pospg(pg));
    endfor
    W = diag(pespg);

    for e = 1:Nelem

      nod = icone(e,:);
      he = xnode(e+1) - xnode(e);
      J = he/2;                                  % Jacobiano (constante, nodo medio centrado)

      % Ke_lm = (k/J) * sum_g w_g dN_l dN_m ;  Me_lm = c*J * sum_g w_g N_l N_m
      Ke = (k/J) * (dN_eval' * W * dN_eval);     % = k/(3*he) [7 -8 1; -8 16 -8; 1 -8 7]
      Me = c*J * (N_eval' * W * N_eval);         % = c*he/30 [4 2 -1; 2 16 2; -1 2 4]

      % Fe_l = J * sum_g w_g N_l(xi_g) G(x(xi_g))
      if isa(G, 'function_handle')
        xg = N_eval * xdof(nod);                 % x(xi_g) = sum_i N_i(xi_g) x_i
        Gg = arrayfun(G, xg);
      else
        Gg = N_eval * G_vec(nod);                % G interpolado en el elemento (constante -> G)
      endif
      Fe = J * (N_eval' * (pespg(:) .* Gg));     % constante: G*he/6*[1;4;1]

      A(nod, nod) += (Ke + Me);
      b(nod) += Fe;

      if esquema ~= 0
        C(nod, nod) += rho*cp*J * (N_eval' * W * N_eval);   % = rho*cp*he/30 [4 2 -1; 2 16 2; -1 2 4]
      endif
    endfor

  endif

  % --- Condiciones de borde (los extremos son nodos de esquina en ambos tipos de elemento)
  fijo = [];    % nodos con Dirichlet
  Tfijo = [];   % valores impuestos en esos nodos

  % Condición de borde izquierdo:
  if cb(1,1) == 3   % Robin
    A(1,1) += cb(1,2);
    b(1) += cb(1,2)*cb(1,3);

  elseif cb(1,1) == 2   % Neumann
    b(1) -= cb(1,2);

  else    % Dirichlet
    A(1,1) = 1;
    A(1, 2:end) = 0;
    b(1) = cb(1,2);
    fijo(end+1) = 1;
    Tfijo(end+1) = cb(1,2);

  endif

  % Condición de borde derecho:
  if cb(2,1) == 3   % Robin
    A(end,end) += cb(2,2);
    b(end) += cb(2,2)*cb(2,3);

  elseif cb(2,1) == 2   % Neumann
    b(end) -= cb(2,2);

  else    % Dirichlet
    A(end,end) = 1;
    A(end, 1:end-1) = 0;
    b(end) = cb(2,2);
    fijo(end+1) = Ndof;
    Tfijo(end+1) = cb(2,2);

  endif

  if esquema == 0    % Estacionario

    T = A\b;
    tvec = [];

  else               % Transiente

    dt    = et(2);
    tol   = et(3);
    maxIt = et(4);
    t0    = et(5);
    tF    = et(6);

    % Esquema temporal -> theta:  (C + theta*dt*A) T^{n+1} = (C - (1-theta)*dt*A) T^n + dt*b
    switch esquema
      case 1
        theta = 0;    % Euler explícito
      case 2
        theta = 1;    % Euler implícito
      otherwise
        error('elemFinitos: esquema temporal %d no implementado (0, 1 o 2).', esquema);
    endswitch

    Tfijo = Tfijo(:);
    libre = setdiff(1:Ndof, fijo);   % nodos incógnita

    % Condición inicial
    if isfield(model, 'T0')
      T0 = model.T0;
      if isa(T0, 'function_handle')
        Tn = arrayfun(T0, xdof);
      elseif length(T0) == Ndof
        Tn = T0(:);
      else
        Tn = T0 * ones(Ndof, 1);
      endif
    else
      Tn = zeros(Ndof, 1);
    endif
    Tn(fijo) = Tfijo;    % compatibilizar con Dirichlet

    % Chequeo de estabilidad (solo esquemas condicionalmente estables, theta < 1/2):
    % dt_crit = 2 / ((1 - 2*theta) * lambda_max),  con  A phi = lambda C phi
    if theta < 0.5
      lambda_max = max(real(eig(C(libre,libre) \ A(libre,libre))));
      dt_crit = 2 / ((1 - 2*theta) * lambda_max);
      if dt > dt_crit
        warning('elemFinitos:inestable', ...
                'dt = %g supera el dt crítico = %g: el esquema explícito es inestable.', dt, dt_crit);
      endif
    endif

    % Matrices constantes del esquema
    Kimp = C + theta*dt*A;            % lado izquierdo (explícito: Kimp = C)
    R    = C - (1 - theta)*dt*A;      % lado derecho   (implícito: R    = C)
    Kll  = Kimp(libre, libre);
    Klf  = Kimp(libre, fijo);

    Npasos = min(maxIt, round((tF - t0)/dt));
    T = zeros(Ndof, Npasos+1);
    tvec = zeros(1, Npasos+1);
    T(:,1) = Tn;
    tvec(1) = t0;

    for n = 1:Npasos
      rhs = R*Tn + dt*b;
      Tnew = Tn;                                   % conserva los valores de Dirichlet
      Tnew(libre) = Kll \ (rhs(libre) - Klf*Tfijo);

      T(:, n+1) = Tnew;
      tvec(n+1) = t0 + n*dt;

      dif = norm(Tnew - Tn, inf);
      Tn = Tnew;

      if dif < tol    % llegó al estado estacionario
        T = T(:, 1:n+1);
        tvec = tvec(1:n+1);
        break;
      endif
    endfor

  endif

endfunction
