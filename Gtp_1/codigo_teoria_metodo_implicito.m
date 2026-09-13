% DifRobinImplicito.m
tic
L = 1; T = 0.5; k = 1;
f = @(x,t)(10*exp(-100*(x-.5).^2).*(t<1.0));  % fuente hasta t=1
a = @(t)(0);
H1 = 3; H2 = 3; uE = @(t)(0);
u0 = @(x)(zeros(size(x)));

N = 30;
h = L/N;

% fijamos deltat y recalculamos lambda correctamente
%deltat = 0.005;
deltat = 0.05;
lambda = k * deltat / h^2; % importante: consistente con deltat

% Condicion inicial
X = linspace(0,L,N+1)';       % nodos fisicos 1..N+1
U = u0(X);                    % U(1:N+1)
% añadimos nodo fantasma U(N+2) inicial usando la formula de Robin
U(N+2) = (H2*uE(0)-H1*U(N+1))*2*h/k + U(N);

ejes = [0 L 0 1];
figure(1); plot(X,U(1:N+1));
title('t = 0');
axis(ejes); pause(0.01)

% Armado de la matriz
%(tamaño N+2 para incluir la fila de la ecuacion Robin)
unos = ones(N+2,1);
columnas = [-lambda*unos (1+2*lambda)*unos -lambda*unos];
matriz = spdiags(columnas, [-1 0 1], N+2, N+2);

% fila Dirichlet en x=0
matriz(1,1:3) = [1 0 0];
% fila Robin: -U_N + (2*h*H1/k)*U_{N+1} + 1*U_{N+2} = 2*h*H2/k * uE
% (coincide con tu construcción original)
matriz(N+2, N:N+2) = [-1 2*h*H1/k 1];

% Resolución temporal
F = zeros(N+2,1);
for t = deltat:deltat:T
    % condiciones en el vector U en tiempo n (lado derecho)
    U(1) = a(t); % Dirichlet en x=0 evaluado en t (estás usando BE -> t^{n+1})
    % NO sobrescribir U(N+2) aquí.
    %La ecuacion Robin está representada en la matriz.
    % Construyo el término fuente en t^{n+1} (Backward Euler)
    F = [0; f(X(2:N+1), t); 0];   % f para nodos 2..N+1, y cero en filas de BCs

    lado_derecho = U + deltat * F;

    % Ahora ajustar el término RHS para la fila de Robin:
    % Recordá que la ecuación de Robin que pusimos en la matriz tiene
    %RHS = 2*h*H2/k*uE(t)
    lado_derecho(N+2) = 2*h*H2/k * uE(t);

    % Resolver sistema implícito
    U = matriz \ lado_derecho;

    % graficar
    figure(1);
    plot(X, U(1:N+1), '*-');
    title(sprintf('Implicito: tiempo fisico = %5.3f', t));
    axis(ejes);
    grid on; grid minor
    pause(0.01);
end
disp('Método Implicito: Tiempo de calculo')
tiempoImplicito = toc
disp('Valor de U al final')
U(end)

