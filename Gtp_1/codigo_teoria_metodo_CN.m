% difusionrobindf1dCN.m
%
% Programita para resolver la ecuacion de Difusion
% con condiciones de Dirichlet en x=0 y Robin en x=L.
% Metodo: Diferencias finitas en espacio y CRANK-NICOLSON en tiempo.
%
% u_t(x,t) - k u_xx(x,t) = f(x,t),    0 < x < L,   t > 0,
% u(0,t) = a(t),     k u'(L,t) + H1 u(L,t) = H2 uE(t),   t > 0,
% u(x,0) = u0(x),    0 < x < L.
% medimos el tiempo de calculo
tic
%% Parametros del problema
L = 1;   T = 0.5;   k = 1;

f = @(x,t)(10*exp(-100*(x-.5).^2)*(t<1.5));
a = @(t)(0);
H1 = 3;   H2 = 3;   uE = @(t)(0);
u0 = @(x)(zeros(size(x)));

%% Parametros del metodo de resolucion
N = 30;
h = L/N;

lambda = 0.25; % el mismo parametro del implicito, solo para fijar deltat
deltat = lambda*h^2/k;

lam = k*deltat/(2*h^2);   % lambda "de Crank-Nicolson" = lambda/2

%% Condicion inicial
X = linspace(0,L,N+1)';
U = u0(X);
% le damos valor a la U en el nodo ficticio
U(N+2) = (H2*uE(0)-H1*U(N+1))*2*h/k + U(N+1);

ejes = [0 L 0 1];
figure(1);   plot(X,U(1:N+1));
title('t = 0');
axis(ejes);   %pause(0.01)

%% Armado de la matriz del lado izquierdo (incognitas en j+1)
unos = ones(N+2,1);
columnas = [-lam*unos (1+2*lam)*unos -lam*unos];

matriz = spdiags(columnas, [-1 0 1], N+2, N+2);
matriz(1,1:3)     = [1 0 0];
matriz(N+2, N:N+2) = [-1  2*h*H1/k  1];

%% Armado de la matriz del lado derecho (conocidos en j)
columnasB = [lam*unos (1-2*lam)*unos lam*unos];
matrizB = spdiags(columnasB, [-1 0 1], N+2, N+2);
% las filas de borde no se promedian: no aportan al lado derecho via matrizB
matrizB(1,:)   = 0;
matrizB(N+2,:) = 0;

%% Resolucion
F = zeros(N+2,1);

for t = deltat:deltat:T
    % t      = t_{j+1}
    % t-deltat = t_j

    % Ensamblado del lado derecho
    F(1)     = a(t);
    F(2:N+1) = (deltat/2)*( f(X(2:N+1),t) + f(X(2:N+1),t-deltat) );
    F(N+2)   = 2*H2*h/k*uE(t);

    lado_derecho = matrizB*U + F;

    % calculamos la nueva U
    U = matriz \ lado_derecho;

    % graficamos la solucion a tiempo t
    figure(1);
    plot(X,U(1:N+1),'*-');
    title(sprintf('Crank-Nicolson: tiempo fisico = %5.3f',t));
    grid on; grid minor
    axis(ejes);   pause(0.01);
end
disp('Método Crank-Nicolson:Tiempo de calculo')
tiempoCN = toc
disp('Valor de U al final')
U(end)

