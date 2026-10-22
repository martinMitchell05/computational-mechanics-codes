function [PHI_vec, Q_vec] = fdm2d_implicit(K,F,xnode,neighb,model,dt)
% Descripción: módulo para resolver el sistema lineal de ecuaciones utilizando esquema
% temporal implícito.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo de la malla.
% * neighb: matriz de vecindad.
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)
% * dt: paso temporal crítico para método explícito.

% Salida:
% * PHI_vec: matriz solución. Cada elemento del vector representa un valor escalar
%   asociado a cada nodo de la malla, y su posición dentro del vector depende de
%   cómo se especificó cada nodo en xnode. Cada columna representa una iteración
%   del esquema temporal (en total nit columnas).
% * Q_vec: matriz de flujo de calor. Para cada nodo se halla un vector bidimensional
%   de flujo de calor, representado por un par (Qx,Qy). Cada par de columnas
%   representa una iteración del esquema temporal (en total 2×nit columnas).
% ----------------------------------------------------------------------

    A = dt / (model.rho * model.cp);

    I = eye(model.nnodes, model.nnodes);

    % Matriz del sistema implícito: A_imp = (I + A * K)
    A_imp = I + A * K;

    % PHI_n := condición inicial del sistema
    PHI_n = model.PHI_n;
    PHI_vec = PHI_n;

    % Flujo de calor inicial en t = 0 ---> parece no ser necesario.
%    [Q_0] = fdm2d_flux(PHI_n, neighb, xnode, model.k);
%    Q_vec = Q_0;
    Q_vec = zeros(model.nnodes,2);


    for n = 1 : model.maxit

        b_imp = PHI_n + A * F;

        % A_imp * PHI = b_imp
        PHI = A_imp \ b_imp;

        % Error relativo entre las últimas dos iteraciones
        err = norm(PHI - PHI_n, 2) / norm(PHI, 2);

        PHI_n = PHI;
        PHI_vec = [PHI_vec, PHI];

        % Cálculo del flujo térmico en el paso actual
        [Q] = fdm2d_flux(PHI, neighb, xnode, model.k);
        Q_vec = [Q_vec, Q];

        % Verificación de tolerancia de corte por estado estacionario
        if err < model.tol
            disp('Método terminado por tolerancia de error.');
            return;
        end
    end

    disp('Método terminado por límite de iteraciones.');

end
