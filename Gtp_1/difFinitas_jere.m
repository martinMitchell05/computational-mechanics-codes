function [T] = difFinitas_jere(xnode, model, cb, et)
  N = length(xnode);

  K = zeros(N,N);
  for i=2:N-1
    dxp = xnode(i+1)-xnode(i);
    dxn = xnode(i)-xnode(i-1);
    a = 2/(dxp*(dxp+dxn));
    b = -2/(dxp*dxn);
    c = 2/(dxn*(dxp+dxn));
    K(i,i-1:i+1) = [model.k*c (model.k*b-model.c) model.k*a];
  endfor

  G = zeros(N,1);
  G(2:N-1) = model.G(2:end-1);

  #Borde izquierdo
  dx_izq = xnode(2)-xnode(1);
  if cb(1,1) == 2 #Neumann
    K(1,1:2) = [(-2*model.k/dx_izq^2 - model.c) 2*model.k/dx_izq^2];
    G(1) = model.G(1) - 2*cb(1,2)/dx_izq;
  elseif cb(1,1) == 3 #Robin
    K(1,1:2) = [(-2*model.k/dx_izq^2 - model.c - 2*cb(1,2)/dx_izq) 2*model.k/dx_izq^2];
    G(1) = model.G(1) + 2*cb(1,2)*cb(1,3)/dx_izq;
  endif

  #Borde derecho
  dx_der = xnode(N)-xnode(N-1);
  if cb(2,1) == 2 #Neumann
    K(N,N-1:N) = [2*model.k/dx_der^2 (-2*model.k/dx_der^2 - model.c)];
    G(N) = model.G(N) - 2*cb(2,2)/dx_der;
  elseif cb(2,1) == 3 #Robin
    K(N,N-1:N) = [2*model.k/dx_der^2 (-2*model.k/dx_der^2 - model.c - 2*cb(2,2)/dx_der)];
    G(N) = model.G(N) + 2*cb(2,2)*cb(2,3)/dx_der;
  endif

  if et(1) == 0 ###  estacionario  ###
    F = -G;
    if cb(1,1) == 1 %Dirichlet
      K(1,1) = 1;
      F(1) = cb(1,2);
    endif
    if cb(2,1) == 1 %Dirichlet
      K(N,N) = 1;
      F(N) = cb(2,2);
    endif
    T = K\F;

  elseif et(1) == 1 ### forward euler ###
    maxIt = et(2);
    tol = et(3);
    alpha = model.k / (model.cp*model.ro);

    dx_min = min(diff(xnode));
    dt = (dx_min^2 / (2*alpha)) * 0.9;

    Tant = zeros(N,1);

    if cb(1,1) == 1 %Dirichlet
      Tant(1) = cb(1,2);
    endif
    if cb(2,1) == 1 %Dirichlet
      Tant(N) = cb(2,2);
    endif

    Thist = [];
    for i = 1:maxIt
      T = ( dt/(model.ro*model.cp) ) * (K*Tant + G) + Tant;

      if cb(1,1) == 1 %Dirichlet
        T(1) = cb(1,2);
      endif
      if cb(2,1) == 1 %Dirichlet
        T(N) = cb(2,2);
      endif

      e = norm(T-Tant)/(norm(Tant));
      if e <= tol
        disp("Se alcanzo la tolerancia");
        break;
      endif
      Thist = [Thist Tant];
      Tant = T;
      if i == maxIt
        disp("Se alcanzo la iteracion maxima");
      endif
    endfor

  elseif et(1) == 2 ### backward euler ###
    maxIt = et(2);
    tol = et(3);
    dt = et(4);
    Tant = zeros(N,1);

    if cb(1,1) == 1 %Dirichlet
      Tant(1) = cb(1,2);
    endif
    if cb(2,1) == 1 %Dirichlet
      Tant(N) = cb(2,2);
    endif

    A = ((model.ro*model.cp)/dt)*eye(N) - K;

    Thist = [];
    for i = 1:maxIt
      b = G + ((model.ro*model.cp)/dt)*Tant;
      T = A \ b;

      if cb(1,1) == 1 %Dirichlet
        T(1) = cb(1,2);
      endif
      if cb(2,1) == 1 %Dirichlet
        T(N) = cb(2,2);
      endif

      e = norm(T-Tant)/(norm(Tant) + eps);
      if e <= tol
        disp("Se alcanzo la tolerancia");
        break;
      endif
      Thist = [Thist Tant];
      Tant = T;

      if i == maxIt
        disp("Se alcanzo la iteracion maxima");
      endif
    endfor

  elseif et(1) == 3 ### crank-nicolson ###
    maxIt = et(2);
    tol = et(3);
    dt = et(4);
    Tant = zeros(N,1);

    if cb(1,1) == 1 %Dirichlet
      Tant(1) = cb(1,2);
    endif
    if cb(2,1) == 1 %Dirichlet
      Tant(N) = cb(2,2);
    endif
    A = ((model.ro*model.cp)/dt)*eye(N) - 0.5 * K;

    Thist = [];
    for i = 1:maxIt
      b = (0.5 * K + ((model.ro*model.cp)/dt)*eye(N)) * Tant + G;
      T = A \ b;

      if cb(1,1) == 1 %Dirichlet
        T(1) = cb(1,2);
      endif
      if cb(2,1) == 1 %Dirichlet
        T(N) = cb(2,2);
      endif

      e = norm(T-Tant)/(norm(Tant) + eps);
      if e <= tol
        disp("Se alcanzo la tolerancia");
        break;
      endif
      Thist = [Thist Tant];
      Tant = T;

      if i == maxIt
        disp("Se alcanzo la iteracion maxima");
      endif
    endfor
  endif

endfunction
