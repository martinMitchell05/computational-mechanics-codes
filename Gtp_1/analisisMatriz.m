function analisisMatriz(A,b)
% caracterizacion de la matriz de un sistema Ax=b
% analisisMatriz.m
% Ejemplo:
% Definir una matriz mal condicionada
%A = [1 2 3; 4 5.0001 6; 7 8 9]; % Casi singular

% Calcular el número de condición en distintas normas
kappa_2 = cond(A);          % Norma 2 (espectral)
kappa_1 = cond(A,1);        % Norma 1 (suma de columnas)
kappa_inf = cond(A,inf);    % Norma infinito (suma de filas)

% Mostrar resultados
fprintf('Número de condición en norma 2: %e\n', kappa_2);
fprintf('Número de condición en norma 1: %e\n', kappa_1);
fprintf('Número de condición en norma infinito: %e\n', kappa_inf);

% Definir un sistema Ax = b y resolverlo
%b = [1; 2; 3]; % Vector de términos independientes
x = A\b; % Resolver sistema

% Calcular el error relativo para evaluar sensibilidad
perturbacion = 1e-6 * randn(size(b)); % Pequeña perturbación en b
b_perturbado = b + perturbacion;
x_perturbado = A\b_perturbado; % Resolver con b perturbado

error_relativo = norm(x_perturbado - x) / norm(x);
error_b = norm(perturbacion) / norm(b);

% Mostrar sensibilidad
fprintf('Error relativo en x debido a perturbación en b: %e\n', error_relativo);
fprintf('Error relativo en b: %e\n', error_b);
fprintf('Amplificación del error por el número de condición: %e\n', error_relativo / error_b);

