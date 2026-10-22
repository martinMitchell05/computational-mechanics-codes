function [K,F] = fdm2d_gen_system(K,F,xnode,neighb,k,c,G)
% Descripción: módulo para ensamblar los términos difusivo, reactivo y fuente de
% todos los nodos de la malla, generando el stencil adecuado dependiendo de si es
% un nodo interior o de frontera.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * xnode: matriz de pares (x,y) representando cada nodo de la malla.
% * neighb: matriz de vecindad.
% * k: conductividad térmica del material. Es un vector que permite representar k(x,y).
% * c: constante de reacción del material. Es un vector que permite representar c(x,y).
% * G: fuente de calor. Es un vector que permite representar G(x,y).

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego del ensamble.
% * F: vector de flujo térmico con modificaciones luego del ensamble.
% ----------------------------------------------------------------------

    M = length(F);
    for p = 1 : M
    
        S = neighb(p, 1);
        E = neighb(p, 2);
        N = neighb(p, 3);
        W = neighb(p, 4);
    
        % calcular distancias entre los nodos vecinos y el nodo p
        if (S ~= -1), dS = abs(xnode(S,2) - xnode(p,2)); end
        if (E ~= -1), dE = abs(xnode(E,1) - xnode(p,1)); end
        if (N ~= -1), dN = abs(xnode(N,2) - xnode(p,2)); end
        if (W ~= -1), dW = abs(xnode(W,1) - xnode(p,1)); end
    
        % coeficientes en eje-x
        if (E == -1)
            ax = 0;               % Nodo ficticio Este
            bx = -2 / (dW * dW);  % Nodo central
            cx = 2 / (dW * dW);   % Nodo interior Oeste
        elseif (W == -1)
            ax = 2 / (dE * dE);   % Nodo interior Este
            bx = -2 / (dE * dE);  % Nodo central
            cx = 0;               % Nodo ficticio Oeste
        else % nodo interior
            ax = 2 / (dE * (dE + dW));
            bx = -2 / (dE * dW);
            cx = 2 / (dW * (dE + dW));
        end
    
        % coeficientes eje-y
        if (N == -1)
            ay = 0;               % Nodo ficticio Norte
            by = -2 / (dS * dS);  % Nodo central
            cy = 2 / (dS * dS);   % Nodo interior Sur
        elseif (S == -1)
            ay = 2 / (dN * dN);   % Nodo interior Norte
            by = -2 / (dN * dN);  % Nodo central
            cy = 0;               % Nodo ficticio Sur
        else % nodo interior
            ay = 2 / (dN * (dN + dS));
            by = -2 / (dN * dS);
            cy = 2 / (dS * (dN + dS));
        end
    
        
        K(p, p) = c(p) - k(p)*bx - k(p)*by;
        F(p) = G(p);
    
        if (S ~= -1), K(p,S) = -k(p)*cy; end
        if (E ~= -1), K(p,E) = -k(p)*ax; end
        if (N ~= -1), K(p,N) = -k(p)*ay; end
        if (W ~= -1), K(p,W) = -k(p)*cx; end

    end

end
