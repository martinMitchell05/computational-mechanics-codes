% Ejercicio 1. a

  L1 = 0;
  L2 = 1;

  T_0 = 10;
  T_3 = 50;

  pcp = 0;
  k = 2;
  c = 0;
  G = 100;

  h = 1/3;

  K = diag(ones(3,1), 1) + diag(-2*ones(4,1)) + diag(ones(3,1), -1)

  K(1,1:2) = [1 0];
  K(end, end-1:end) = [0 1];

  disp(K)

  f = [T_0; -(G/k)*(h^2); -(G/k)*(h^2); T_3]

  printf("Solución al sistema por diferencia finita:\n")
  T = K\f;
  disp(T)


  % Este ejemplo se pudo despejar fácil analíticamente las incógnitas.

