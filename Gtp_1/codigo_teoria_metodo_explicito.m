% difusionrobindf1dexpl.m
%
% Programa para resolver la ecuacion de Difusion
% con condiciones de Dirichlet en x=0 y Robin o Mixtas en x=L.
% Metodo: Diferencias finitas en espacio y Euler EXPLICITO en tiempo.
%
% u_t(x,t) - k u_xx(x,t) = f(x,t), 0 < x < L, t > 0
% u(0,t) = a(t), k u’(L,t) + H1 u(L,t) = H2 uE(t), t > 0
% u(x,0) = u0(x), 0 < x < L.
% medimos el tiempo de calculo

tic
% % Parametros del problema
L = 1; T = 0.5; k = 1; % difusividad
f = @(x,t)(10*exp(-100*(x-.5).^2)*(t<0.5));
a = @(t)(0);

% Suponemos iguales los valores del parametro H de Robin
H1 = 3; H2 = 3; uE = @(t)(0); %Temperatura en el exterior
u0 = @(x)(zeros(size(x))); % valor inicial
% % Parametros del metodo de resolucion
N = 30;
h = L/N;
lambda = 0.25; % la mitad de lo necesario para tener estabilidad
deltat = lambda*h^2/k;


% % Condicion inicial
X = linspace(0,L,N+1)';
U = u0(X);
% le damos valor a la U en el nodo ficticio (t=0)
U(N+2) = (H2*uE(0)-H1*U(N+1))*2*h/k + U(N);
ejes = [0 L 0 1];
figure(1); plot(X,U(1:N+1));
title('t = 0')
axis(ejes); pause(0.01)


% % Resolucion
for t=deltat:deltat:T
  % calculamos U en los nodos "interiores"
  U(2:N+1) = deltat*f(X(2:N+1), t) ...
  + lambda*U(1:N) + (1-2*lambda)*U(2:N+1) + lambda*U(3:N+2);
  % le damos valor a U en el extremo izquierdo
  U(1) = a(t);
  % le damos valor a la U en el nodo ficticio
  U(N+2) = (H2*uE(t)-H1*U(N+1))*2*h/k + U(N);
  % graficamos la solucion a tiempo t
  figure(1);
  plot(X,U(1:N+1),'*-');
  title(sprintf('Explicito: tiempo fisico = %5.3f',t));
  axis(ejes);
  grid on; grid minor
  pause(0.01);
end


disp('Método Explicito: Tiempo de calculo')
tiempoExplicito = toc
disp('valor de U al final')
U(end)


