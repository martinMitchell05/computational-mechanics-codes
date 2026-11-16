function [K,F] = fdm2d_robin(K,F,xnode,neighb,ROB)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos 
% pertenecientes a fronteras de tipo Robin.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de 
%   cada nodo de la malla.
% * neighb: matriz de vecindad.
% * ROB: matriz con la información sobre la frontera de tipo Robin.
%   - Columna 1: índice del nodo donde se aplica la condición de borde.
%   - Columna 2: valor de coeficiente de calor (h)
%   - Columna 3: valor de temperatura de referencia (phi_inf).
%   - Columna 4: dirección y sentido del flujo:
%     (1) Flujo en dirección eje-y, sentido negativo (S – South – Sur)
%     (2) Flujo en dirección eje-x, sentido positivo (E – East – Este)
%     (3) Flujo en dirección eje-y, sentido positivo (N – North – Norte)
%     (4) Flujo en dirección eje-x, sentido negativo (W – West – Oeste)

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego de 
%   aplicar la condición de borde.
% * F: vector de flujo térmico con modificaciones luego de aplicar la condición 
%   de borde.
% ----------------------------------------------------------------------

    for i = 1 : size(ROB, 1)

        nodo = ROB(i, 1);
        S = neighb(nodo, 1);
        E = neighb(nodo, 2);
        N = neighb(nodo, 3);
        W = neighb(nodo, 4);

        h = ROB(i, 2);
        phi_inf = ROB(i, 3);
        
        switch (ROB(i, 4))
            case 1
                dy = abs(xnode(N, 2) - xnode(nodo, 2)); % para calcular la distancia en "y" usar la columna 2
                F(nodo) += 2*h*phi_inf/dy;
                K(nodo, nodo) += 2*h/dy;     % A lo que habia en el nodo '(i,j)' se le adiciona lo que sale de condicion de Robin 
                
            case 2
                dx = abs(xnode(W, 1) - xnode(nodo, 1));
                F(nodo) += 2*h*phi_inf/dx;
                K(nodo, nodo) += 2*h/dx;
                
            case 3
                dy = abs(xnode(S, 2) - xnode(nodo, 2));
                F(nodo) += 2*h*phi_inf/dy;
                K(nodo, nodo) += 2*h/dy;
                
            case 4
                dx = abs(xnode(E, 1) - xnode(nodo, 1));
                F(nodo) += 2*h*phi_inf/dx;
                K(nodo, nodo) += 2*h/dx;
        end
    end

end