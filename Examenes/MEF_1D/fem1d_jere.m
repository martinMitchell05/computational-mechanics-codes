% =========================================================================
% Código 1D de Elementos Finitos (Galerkin) en GNU Octave
% Soporta elementos lineales (2 nodos) y cuadráticos (3 nodos)
% Soporta término fuente G(x) variable
% Condiciones de borde: Dirichlet, Neumann y Robin
% Simulación estacionaria o transitoria (Euler explícito/implícito)
% =========================================================================

clear all; close all; clc;

% -------------------------------------------------------------------------
% 1. CONFIGURACIÓN DEL TIPO DE ELEMENTO Y LA MALLA
% -------------------------------------------------------------------------
tipo_elemento = 'cuadratico'; % Opciones: 'lineal', 'cuadratico'

% Para elementos cuadráticos, los nodos intermedios deben estar en el punto
% medio geométrico de cada elemento para que las matrices elementales
% cuadráticas usadas en este programa sean válidas.
X = [0.0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0];
Nnodos = length(X);

% Cálculo de la cantidad de elementos según el tipo elegido
if strcmp(tipo_elemento, 'lineal')
    Nelem = Nnodos - 1;
elseif strcmp(tipo_elemento, 'cuadratico')
    if mod(Nnodos, 2) == 0
        error('Para elementos cuadráticos, el número total de nodos de X debe ser impar.');
    end
    Nelem = (Nnodos - 1) / 2;
else
    error('Tipo de elemento no soportado. Elija "lineal" o "cuadratico".');
end

% Validar que las coordenadas estén estrictamente ordenadas
if any(diff(X) <= 0)
    error('Las coordenadas de X deben ser estrictamente crecientes.');
end

% -------------------------------------------------------------------------
% 2. PARÁMETROS DEL PROBLEMA FÍSICO
% -------------------------------------------------------------------------
k     = 1.0;   % Coeficiente de difusión/conductividad
c     = 0.0;   % Coeficiente de reacción
rhoCp = 2.0;   % Capacidad temporal (rho * Cp), usada en transitorio

% Término fuente G(x)
G_func = @(x) 50.0;

% -------------------------------------------------------------------------
% 3. CONFIGURACIÓN DE CONDICIONES DE BORDE
% -------------------------------------------------------------------------
% Tipos soportados: 'dirichlet', 'neumann', 'robin'
%
% Convención de signos:
% - Neumann: cb_valor1 representa k * d(phi)/dn, con n normal exterior.
% - Robin: -k*d(phi)/dn = h*(phi - phi_infinito).

cb_izq_tipo   = 'robin';
cb_izq_valor1 = 10;   % h (Robin), k*d(phi)/dn (Neumann) o phi (Dirichlet)
cb_izq_valor2 = 5;    % phi_infinito (Robin)

cb_der_tipo   = 'dirichlet';
cb_der_valor1 = 25;   % Valor fijo de Dirichlet
cb_der_valor2 = 0;    % No aplica

% -------------------------------------------------------------------------
% 4. CONFIGURACIÓN DE LA SIMULACIÓN
% -------------------------------------------------------------------------
tipo_simulacion = 'transitorio'; % 'estacionario' o 'transitorio'
metodo_temporal = 'implicito';   % 'explicito' o 'implicito'
t_final   = 2.0;
dt        = 0.05;
T_inicial = 0.0;

% -------------------------------------------------------------------------
% 5. INICIALIZACIÓN Y ENSAMBLE DE MATRICES
% -------------------------------------------------------------------------
K = zeros(Nnodos, Nnodos); % Difusión
M = zeros(Nnodos, Nnodos); % Reacción
C = zeros(Nnodos, Nnodos); % Capacidad temporal (matriz de masa)
F = zeros(Nnodos, 1);      % Fuente

for e = 1:Nelem

    if strcmp(tipo_elemento, 'lineal')
        % --- FORMULACIÓN LINEAL (2 NODOS) ---
        n1 = e;
        n2 = e + 1;
        he = X(n2) - X(n1);
        indices = [n1, n2];

        Ke = (k / he) * [1, -1; -1, 1];
        Me = (c * he / 6) * [2, 1; 1, 2];
        Ce = (rhoCp * he / 6) * [2, 1; 1, 2];

        % Fuente integrada con interpolación lineal de los valores nodales
        G1 = G_func(X(n1));
        G2 = G_func(X(n2));
        Fe = (he / 6) * [2, 1; 1, 2] * [G1; G2];

    elseif strcmp(tipo_elemento, 'cuadratico')
        % --- FORMULACIÓN CUADRÁTICA (3 NODOS) ---
        n1 = 2*e - 1;
        n2 = 2*e;
        n3 = 2*e + 1;
        he = X(n3) - X(n1); % Longitud total del elemento
        indices = [n1, n2, n3];

        % Estas matrices suponen que X(n2) = (X(n1) + X(n3))/2.
        % La malla X definida arriba cumple esta condición en cada elemento.
        Ke = (k / (3 * he)) * [ 7, -8,  1;
                               -8, 16, -8;
                                1, -8,  7 ];

        Me = (c * he / 30) * [ 4,  2, -1;
                               2, 16,  2;
                              -1,  2,  4 ];

        Ce = (rhoCp * he / 30) * [ 4,  2, -1;
                                   2, 16,  2;
                                  -1,  2,  4 ];

        % Fuente interpolada cuadráticamente sobre los tres nodos
        G1 = G_func(X(n1));
        G2 = G_func(X(n2));
        G3 = G_func(X(n3));
        Fe = (he / 30) * [ 4,  2, -1;
                           2, 16,  2;
                          -1,  2,  4 ] * [G1; G2; G3];
    end

    % Ensamble global
    K(indices, indices) = K(indices, indices) + Ke;
    M(indices, indices) = M(indices, indices) + Me;
    C(indices, indices) = C(indices, indices) + Ce;
    F(indices) = F(indices) + Fe;
end

A_global = K + M;
b_global = F;

% -------------------------------------------------------------------------
% 6. APLICACIÓN DE CONDICIONES DE NEUMANN Y ROBIN
% -------------------------------------------------------------------------
H = zeros(Nnodos, Nnodos);
Q = zeros(Nnodos, 1);

if strcmp(cb_izq_tipo, 'neumann')
    Q(1) = Q(1) + cb_izq_valor1;
elseif strcmp(cb_izq_tipo, 'robin')
    H(1, 1) = H(1, 1) + cb_izq_valor1;
    Q(1)    = Q(1) + cb_izq_valor1 * cb_izq_valor2;
end

if strcmp(cb_der_tipo, 'neumann')
    Q(Nnodos) = Q(Nnodos) + cb_der_valor1;
elseif strcmp(cb_der_tipo, 'robin')
    H(Nnodos, Nnodos) = H(Nnodos, Nnodos) + cb_der_valor1;
    Q(Nnodos)         = Q(Nnodos) + cb_der_valor1 * cb_der_valor2;
end

A_global = A_global + H;
b_global = b_global + Q;

% -------------------------------------------------------------------------
% 6.1. CONTROL DEL PASO TEMPORAL PARA EULER EXPLÍCITO
% -------------------------------------------------------------------------
% La fórmula dt <= dx^2/(2*alpha) corresponde a diferencias finitas 1D
% explícitas con malla uniforme. Aquí se usa dx_min como cota conservadora
% orientativa porque la malla no es uniforme y el método espacial es FEM.
%
% Para FEM también se calcula el límite espectral de Euler hacia adelante:
% dt <= 2/lambda_max, donde A*v = lambda*C*v sobre los grados de libertad
% que no están fijados por condiciones Dirichlet.
%
% Se toma el menor de ambos límites con un margen de seguridad del 10 %.
if strcmp(tipo_simulacion, 'transitorio') && strcmp(metodo_temporal, 'explicito')
    if k <= 0 || rhoCp <= 0
        error('Para estimar la estabilidad explícita, k y rhoCp deben ser positivos.');
    end

    alpha = k / rhoCp;
    dx_min = min(diff(X));
    dt_limite_dif_fin = dx_min^2 / (2 * alpha);

    grados_libres = 1:Nnodos;
    if strcmp(cb_izq_tipo, 'dirichlet')
        grados_libres(grados_libres == 1) = [];
    end
    if strcmp(cb_der_tipo, 'dirichlet')
        grados_libres(grados_libres == Nnodos) = [];
    end

    if isempty(grados_libres)
        error('No hay grados de libertad libres para calcular el límite de estabilidad FEM.');
    end

    A_libre = A_global(grados_libres, grados_libres);
    C_libre = C(grados_libres, grados_libres);
    autovalores = eig(A_libre, C_libre);
    lambda_max = max(real(autovalores));

    if isempty(lambda_max) || ~isfinite(lambda_max) || lambda_max <= 0
        error('No se pudo calcular un límite de estabilidad explícita válido a partir de las matrices FEM.');
    end

    dt_limite_fem = 2 / lambda_max;
    factor_seguridad = 0.90;
    dt_recomendado = factor_seguridad * min(dt_limite_dif_fin, dt_limite_fem);

    fprintf('\n--- Control de estabilidad de Euler explicito ---\n');
    fprintf('Difusividad termica alpha = %.8g\n', alpha);
    fprintf('Espaciamiento minimo dx_min = %.8g\n', dx_min);
    fprintf('Limite orientativo por formula de diferencias finitas = %.8g s\n', dt_limite_dif_fin);
    fprintf('Limite espectral calculado para FEM = %.8g s\n', dt_limite_fem);
    fprintf('Paso maximo recomendado (margen del 10%%) = %.8g s\n', dt_recomendado);

    if dt > dt_recomendado
        fprintf('El dt solicitado (%.8g s) supera el recomendado; se ajusta a %.8g s.\n', dt, dt_recomendado);
        dt = dt_recomendado;
    else
        fprintf('El dt solicitado (%.8g s) cumple el limite recomendado.\n', dt);
    end
    fprintf('------------------------------------------------\n\n');
end

% -------------------------------------------------------------------------
% 7. RESOLUCIÓN SEGÚN EL TIPO DE SIMULACIÓN
% -------------------------------------------------------------------------
if strcmp(tipo_simulacion, 'estacionario')
    A_eff = A_global;
    b_eff = b_global;

    if strcmp(cb_izq_tipo, 'dirichlet')
        A_eff(1, :) = 0;
        A_eff(1, 1) = 1;
        b_eff(1) = cb_izq_valor1;
    end
    if strcmp(cb_der_tipo, 'dirichlet')
        A_eff(Nnodos, :) = 0;
        A_eff(Nnodos, Nnodos) = 1;
        b_eff(Nnodos) = cb_der_valor1;
    end

    T = A_eff \ b_eff;

    figure('Name', 'FEM 1D - Estacionario');
    plot(X, T, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    title(sprintf('Solución Estacionaria (Elementos: %s)', tipo_elemento));
    xlabel('Posición x'); ylabel('Variable \phi');
    grid on;

elseif strcmp(tipo_simulacion, 'transitorio')
    pasos = ceil(t_final / dt); % Se mantiene el comportamiento original solicitado
    T = ones(Nnodos, 1) * T_inicial;

    % Imponer Dirichlet en la condición inicial, antes del primer paso temporal
    if strcmp(cb_izq_tipo, 'dirichlet')
        T(1) = cb_izq_valor1;
    end
    if strcmp(cb_der_tipo, 'dirichlet')
        T(Nnodos) = cb_der_valor1;
    end

    if strcmp(metodo_temporal, 'explicito')
        Mat_LHS = C;
        Mat_RHS = C - dt * A_global;
    elseif strcmp(metodo_temporal, 'implicito')
        Mat_LHS = C + dt * A_global;
        Mat_RHS = C;
    else
        error('Método temporal no reconocido. Use "explicito" o "implicito".');
    end

    % Aplicar Dirichlet en el sistema temporal
    if strcmp(cb_izq_tipo, 'dirichlet')
        Mat_LHS(1, :) = 0;
        Mat_LHS(1, 1) = 1;
    end
    if strcmp(cb_der_tipo, 'dirichlet')
        Mat_LHS(Nnodos, :) = 0;
        Mat_LHS(Nnodos, Nnodos) = 1;
    end

    figure('Name', 'FEM 1D - Transitorio');
    hold on; grid on;
    xlabel('Posición x'); ylabel('Variable \phi');
    title(sprintf('Evolución Transitoria (%s - %s)', metodo_temporal, tipo_elemento));

    for n = 1:pasos
        b_eff = Mat_RHS * T + dt * b_global;

        if strcmp(cb_izq_tipo, 'dirichlet')
            b_eff(1) = cb_izq_valor1;
        end
        if strcmp(cb_der_tipo, 'dirichlet')
            b_eff(Nnodos) = cb_der_valor1;
        end

        T_new = Mat_LHS \ b_eff;
        T = T_new;

        if mod(n, max(1, ceil(pasos/5))) == 0 || n == 1
            plot(X, T, '-o', 'DisplayName', sprintf('t = %.2f s', n*dt));
        end
    end

    legend('Location', 'northwest');
    hold off;
else
    error('Tipo de simulación no reconocido. Use "estacionario" o "transitorio".');
end
