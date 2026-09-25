function [Q] = fdm2d_flux(PHI,neighb,xnode,k)
% Descripción: módulo calcular el flujo de calor en todo el dominio. Se aplica la
% Ley de Fourier y se evalúa como fluye el calor en todos los puntos (nodos) del dominio.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * PHI: vector solución. Cada elemento del vector representa un valor escalar asociado
% a cada nodo de la malla, y su posición dentro del vector depende de cómo se especificó
% cada nodo en xnode.

% Salida:
% * Q: vector de flujo de calor. Para cada nodo se halla un vector bidimensional de 
%   flujo de calor, representado por un par (Qx,Qy)
% ----------------------------------------------------------------------

    Q = [];
end

