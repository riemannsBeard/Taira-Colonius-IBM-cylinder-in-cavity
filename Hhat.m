function [H, beta] = Hhat(grid, ib, Nx, Ny)

    ds.x = ib.dxi';
    ds.y = ib.deta';
    
    beta = hypot(ds.x, ds.y);
    
    %% u-velocity
    r.x = ib.xi' - reshape(grid.xu, [], 1);
    r.y = ib.eta' - reshape(grid.yu, [], 1);
    
    dr.x = repmat(grid.dX, [Ny-1,1]);
    dr.y = repmat(grid.dYu', [Nx,1]);

    H.u = beta.*delta(r.x, dr.x).*delta(r.y, dr.y);

    %% v-velocity
    r.x = ib.xi' - reshape(grid.xv, [], 1);
    r.y = ib.eta' - reshape(grid.yv, [], 1);
    
    dr.x = repmat(grid.dXv, [1,Ny])';
    dr.y = repmat(grid.dY', [1,Nx-1])';
    
    H.v = beta.*delta(r.x, dr.x).*delta(r.y, dr.y);

end

