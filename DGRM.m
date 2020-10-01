function [ D, G, R, M ] = DGRM( grid, Nx, Ny )
% Funcion que devuelve los operadores matriciales G, D, R y M

    % Vectores unitarios auxiliares
    ex = ones(Nx+1, 1);
    ey = ones(Ny+1, 1);

    %% OPERADOR GRADIENTE
    % Se aplica a la presion
        
    % Operador derivada para una fila y una columna
    g.x = spdiags([-ex, ex], [0, 1], Nx-1, Nx);
    g.y = spdiags([-ey, ey], [0, 1], Ny-1, Ny);

    G.x = kron(g.x, speye(Ny));
    G.y = kron(speye(Nx), g.y);

    G.G = [G.x; G.y];
           
    %% OPERADOR DIVERGENCIA
    % Se aplica a las velocidades
    
    D.D = -G.G';
    
    % Condiciones de contorno

    d.uW = spdiags(-ex, 0, Nx, 1);
    D.uW = kron(d.uW, speye(Ny));
    
    d.uE = spdiags(ex, -Nx+1, Nx, 1);
    D.uE = kron(d.uE, speye(Ny));
    
    d.vS = spdiags(-ey, 0, Ny, 1);
    D.vS = kron(speye(Nx), d.vS);
    
    d.vN = spdiags(ey, -Ny+1, Ny, 1);
    D.vN = kron(speye(Nx), d.vN);
        
    %% MATRIZ DE FLUJO
    
    dyj = spdiags(grid.dY, 0, Ny, Ny);
    dxi = spdiags(grid.dX, 0, Nx, Nx);
    
%     R.u = kron(dyj, speye(Nx-1));
%     R.v = kron(speye(Ny-1), dxi);

    R.u = kron(speye(Nx-1), dyj);
    R.v = kron(dxi, speye(Ny-1));

    R.R = blkdiag(R.u, R.v);
           
    %% MATRIZ DE MASA
    
    ix = [0.75; ones(Nx-2,1); 0.75];
    iy = [0.75; ones(Ny-2,1); 0.75];
    
    Ix = spdiags(ix, 0, Nx, Nx);
    Iy = spdiags(iy, 0, Ny, Ny);
    
    Dxp = spdiags(grid.dXp, 0, Nx-1, Nx-1);
    Dyp = spdiags(grid.dYp, 0, Ny-1, Ny-1);

    % Mhat
    Mhat.u = kron(Dxp, Iy);
    Mhat.v = kron(Ix, Dyp);

    M.hat = blkdiag(Mhat.u, Mhat.v);
    
    % M
    M.M = R.R\M.hat;
    
end

