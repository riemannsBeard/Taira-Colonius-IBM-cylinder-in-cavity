clc
clear all
close all

%% Datos

Re = 100;
Nx = 256; % Celdillas en x
Ny = 256; % Celdillas en y
Lx = 1;
Ly = 1;
tf = 4.5;
CFL = 0.5;
Ulid = 1;
Nk = 64;

%% Staggered Grid Generation

alpha.x = 0.;
alpha.y = 0.;

[grid, u, v, p] = gridGeneration(Lx, Ly, Nx, Ny, alpha);
dt = CFL*min(grid.cellMin^2*Re, grid.cellMin);

%% Immersed Boundary Grid Generation

R = 0.25;
ib = IBMcylinder(R, Nk);
ib.xi = ib.xi + 0.5*Lx;
ib.eta = ib.eta + 0.5*Ly;

%% Conventional Boundary Conditions

bc.uS = zeros(1,Nx-1);
bc.uN = ones(1,Nx-1);
bc.uE = zeros(Ny,1);
bc.uW = zeros(Ny,1);

bc.vS = zeros(1,Nx);
bc.vN = zeros(1,Nx);
bc.vE = zeros(Ny-1,1);
bc.vW = zeros(Ny-1,1);

%% Immersed Boundary Conditions

ib.Utheta = pi/8;

ib.u = ib.Utheta*sin(ib.theta);
ib.v = -ib.Utheta*cos(ib.theta);

% figure,
% plot(grid.xc(1,1), grid.yc(1,1), 'kx', grid.xu(1,1), grid.yu(1,1), 'g>',...
%     grid.xv(1,1), grid.yv(1,1), 'r^')
% h = legend('$p$', '$u$', '$v$');
% set(h, 'interpreter', 'latex', 'fontsize', 14)
% %hold on
% plot(grid.x, grid.y, 'k.', grid.x, grid.y, 'k-', grid.y, grid.x, 'k-',...
%     grid.xc, grid.yc, 'kx', grid.xv, grid.yv, 'r^', grid.xu, grid.yu, 'g>')
% set(gca, 'fontname', 'times', 'fontsize', 12);


%% Operators

[D, G, R, M] = DGRM(grid, Nx, Ny);

Lhat = laplacianHat(Nx, Ny, grid);
L.L = M.hat*Lhat.L/R.R;

M.M = M.hat/R.R;
M.inv = inv(M.M);

Ahat = sparse(speye(size(Lhat.L))/dt - 0.5*Lhat.L/Re);
A = M.hat*Ahat/R.R;
dA = decomposition(A);

BN = dt*speye(size(M.M))/M.M + (0.5/Re)*dt*dt*(M.inv*L.L)*M.inv +...
    ((0.5/Re)^2)*(dt^3)*((M.inv*L.L)^2)*M.inv;

%% IBM stuff

% Regularization
[Hhat_, beta] = Hhat(grid, ib, Nx, Ny);
Hhat_ = blkdiag(Hhat_.u, Hhat_.v);
H = sparse(M.M*Hhat_);

% Interpolation
[Ehat_, alpha] = Ehat(grid, ib, Nx, Ny);
Ehat_ = blkdiag(Ehat_.u, Ehat_.v);
E = sparse(Ehat_/R.R);

EH = E*H;

%% Left-Hand Side term

Q = [G.G, E'];

LHS = sparse(Q'*BN*Q);
dLHS = decomposition(LHS);
 
%% Simulation

u = reshape(u, [], 1);
v = reshape(v, [], 1);

uOld = u;
vOld = v;

t = 0;
k = 0;

tic
while t<= tf
    
    u = reshape(u, [], 1);
    v = reshape(v, [], 1);

    % Advective terms
    [NhatOld, ~, ~] = convectionHat(grid, uOld, vOld, Nx, Ny, bc);
    [Nhat, ua, va] = convectionHat(grid, u, v, Nx, Ny, bc);
    
    rnHat = explicitTerms(Lhat, Re, dt, Nhat, NhatOld, u, v);  
    rn = M.hat*rnHat;
        
    %% 1. Solve for intermediate velocity
       
    % BC's due to Laplacian
    bc1hat.u = Lhat.ux0*bc.uW + Lhat.ux1*bc.uE + Lhat.uy1*bc.uN' + ...
        Lhat.uy0*bc.uS';
    bc1hat.v = Lhat.vx0*bc.vW + Lhat.vx1*bc.vE + Lhat.vy1*bc.vN' + ...
        Lhat.vy0*bc.vS';

    bc1 = M.hat*[bc1hat.u; bc1hat.v]/Re;
 
    r1 = rn + bc1;
    
    % Flux calculation
    
    q = dA\r1;    
    qu = q(1:Ny*(Nx-1));
    qv = q(Ny*(Nx-1)+1:end);
       
  
    %% 2. Solve the Poisson Equation
    
    bc2 = D.uW*(bc.uW.*grid.dY) + D.uE*(bc.uE.*grid.dY) + ...
        D.vS*(bc.vS'.*grid.dX) + D.vN*(bc.vN'.*grid.dX);
    
    r2 = [bc2; ib.u; ib.v];

    RHS = Q'*q - r2;
    
    lambda = dLHS\RHS;
    
    % Pressure
    phi = lambda(1:end-2*Nk);
    
    % Forces
    fTilda.x = lambda(end-2*Nk+1:end-Nk);
    fTilda.y = lambda(end-Nk+1:end);
    fTilda.f = [fTilda.x; fTilda.y];
    
    f.f = -(EH)\E*E'*fTilda.f;
    

%     Ehat_ = Ehat(grid, ib, Nx, Ny);
%     Ehat_ = blkdiag(Ehat_.u, Ehat_.v);
%     
%     E = Ehat_/R.R;
%     
%     Q = [G.G, E'];
    
%     BN = dt*speye(size(M.M))/M.M + (0.5/Re)*dt*dt*(M.inv*L.L)*M.inv +...
%         ((0.5/Re)^2)*(dt^3)*((M.inv*L.L)^2)*M.inv;
            
    % LHS = sparse(G.G'*BN*G.G);
    % dLHS = decomposition(LHS);
    
    %LHS = sparse(Q'*BN*Q);
        

    %lambda = LHS\RHS;
    
    %% 3. Projection step
    
    q = q - BN*Q*lambda;
       
    vel = R.R\q;
    
%     ib.xi = ib.xi + ib.u*dt;
%     ib.eta = ib.eta + ib.v*dt;
    
    t = t + dt;
    k = k + 1;

    % Residuals
    epsU(k) = max(abs(u - vel(1:Ny*(Nx-1))));
    epsV(k) = max(abs(v - vel(Ny*(Nx-1)+1:end)));

    % Separation of velocity components
    u = vel(1:Ny*(Nx-1));
    v = vel(Ny*(Nx-1)+1:end);
    
    % Forces storing
    f.x(k) = sum(f.f(1:Nk));
    f.y(k) = sum(f.f(Nk+1:end));
    F(k) = hypot(f.x(k), f.y(k));

    
    % Information
    fprintf(['t = ' num2str(t) '. Elapsed time: ' num2str(toc) 's \n']);
    fprintf(['Residuals u = ' num2str(epsU(k)) '.\n'...
        'Residuals v = ' num2str(epsV(k)) '\n \n']);

%     if (mod(k, 0.25) == 0)
%         u = reshape(u, Ny, Nx-1);
% 
%         F(j) = getframe(u);
%         writeVideo, movie
%     
%     end
        
%     if (mod(k, 50) == 0)
%         u = reshape(u, Ny, Nx-1);
%         figure(2),
%         pcolor(grid.x, grid.y, sqrt(ua.^2 + va.^2)), shading interp
%         colorbar, colormap jet
%         title(['t = ' num2str(t)])
%         drawnow
%     end
    

end
toc

%% Plots
close all

u = reshape(u, Ny, Nx-1);
v = reshape(v, Ny-1, Nx);
phi = reshape(phi, Ny, Nx);
qu = reshape(qu, Ny, Nx-1);
qv = reshape(qv, Ny-1, Nx);

% Contours
figure(2),
pcolor(grid.x, grid.y, hypot(ua,va)), shading interp, hold on
colorbar, colormap jet
%fill(ib.xi, ib.eta, 'w')
plot(ib.xi, ib.eta, 'w.')

title(['$u$; $t = $' num2str(t)], 'interpreter', 'latex',...
    'fontsize', 20)
pbaspect([Lx Ly 1])

figure(3),
s = 2;
quiver(grid.x(1:s:end,1:s:end), grid.y(1:s:end,1:s:end),...
    ua(1:s:end,1:s:end), va(1:s:end,1:s:end), 8), hold on,
fill(ib.xi, ib.eta, 'w')
xlim([0 Lx]), ylim([0 Ly])
pbaspect([Lx Ly 1])

figure(4),
xstart = Lx*rand(Nx,1); 
ystart = Ly*rand(Ny,1); 
streamline(grid.x, grid.y, ua, va, xstart, ystart, [0.03, 3000]);
hold on,
fill(ib.xi, ib.eta, 'w')
xlim([0 Lx]), ylim([0 Ly])
pbaspect([Lx Ly 1])

figure(5),
contourf(grid.xp, grid.yp, phi), hold on, %shading interp
fill(ib.xi, ib.eta, 'w')
colorbar, %colormap jet
title(['$\phi$; $t = $' num2str(t)], 'interpreter', 'latex',...
    'fontsize', 20)
pbaspect([Lx Ly 1])

% Residuals
figure(6),
plot(1:k, epsU, 1:k, epsV)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
h = legend('$u$', '$v$');
set(h, 'interpreter', 'latex', 'fontsize', 16)
xlabel('$N$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\xi$', 'interpreter', 'latex', 'fontsize', 16)
title('Residuals')

% Forces
figure(7),
plot(dt*(1:k), 2*f.x, dt*(1:k), 2*f.y, dt*(1:k), 2*F)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
h = legend('$C_x$', '$C_y$', '$C_F$', 'Location', 'Best');
set(h, 'interpreter', 'latex', 'fontsize', 16)
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 16)
title('Forces')

%% Validation
close all

uData = load('validation/uData_Re1000');
vData = load('validation/vData_Re1000');

% Interpolation
yq = grid.Y; %linspace(0, 1, 64);
xq = yq*0 + 0.5;

uInt = interp2(grid.x, grid.y, ua, xq, yq);
vInt = interp2(grid.x, grid.y, va, yq, xq);

% Comparison

figure(4),
plot(uData(:,1), uData(:,2), 'rs', 'markersize', 8), hold on,
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
xlabel('$y$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$u$', 'interpreter', 'latex', 'fontsize', 16)
plot(uInt, yq)
hold off

figure(5),
plot(vData(:,1), vData(:,2), 'rs', 'markersize', 8), hold on,
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$v$', 'interpreter', 'latex', 'fontsize', 16)
plot(yq, vInt)

