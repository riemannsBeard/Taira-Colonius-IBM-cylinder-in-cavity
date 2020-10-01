function [E, alpha] = Ehat(grid, ib, Nx, Ny)

    %% u-velocity
    r.x = reshape(grid.xu, 1, []) - ib.xi;
    r.y = reshape(grid.yu, 1, []) - ib.eta;
    
    dr.x = repmat(grid.dX', [1,Ny-1]);
    dr.y = repmat(grid.dYu, [1,Nx]);
    
    alpha.u = dr.x.*dr.y;

    E.u = alpha.u.*delta(r.x, dr.x).*delta(r.y, dr.y);

    %% v-velocity
    r.x = reshape(grid.xv, 1, []) - ib.xi;
    r.y = reshape(grid.yv, 1, []) - ib.eta;
    
    dr.x = repmat(grid.dXv', [Ny,1])';
    dr.y = repmat(grid.dY, [Nx-1,1])';
    
    alpha.v = dr.x.*dr.y;

    E.v = alpha.v.*delta(r.x, dr.x).*delta(r.y, dr.y);

end

