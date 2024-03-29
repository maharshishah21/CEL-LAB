% Define heat plate dimensions and properties
Nx = 25; % Number of grid points in x-direction
Ny = 25; % Number of grid points in y-direction
Lx = 1; % Length of the plate in x-direction (m)
Ly = 1; % Length of the plate in y-direction (m)
Lt = 75;% no of time steps
dx = Lx / (Nx - 1); % Grid spacing in x-direction
dy = Ly / (Ny - 1); % Grid spacing in y-direction
dt = 0.001; %time step size
alpha = 0.1; % Thermal diffusivity (m^2/s)
T_const = 50;

% Set initial temperature (hot center, cool edges)
T = 25*ones(Ny, Nx);
T(floor(Nx/3):ceil(2*Nx/3), floor(Ny/3):ceil(2*Ny/3)) = T_const;
              
[X, Y] = meshgrid(linspace(0, L, Nx), linspace(0, W, Ny));
figure;
contourf(X, Y, T', 250, 'LineStyle', 'none');
colorbar;
title('Initial Temperature Field');
xlabel('x');
ylabel('y');

for t = 1:Lt
    T_old = T;
    for i = 2:Nx-1
        for j = 2:Ny-1
T(i,j) = T_old(i,j) + alpha * dt * ((T_old(i+1,j) - 2*T_old(i,j) + T_old(i-1,j)) / dx^2 + ... 
    (T_old(i,j+1) - 2*T_old(i,j) + T_old(i,j-1)) / dy^2);        
        end
    end
% Plot temperature field at selected time steps
    if mod(t, 10) == 0
        figure;
        contourf(X, Y, T', 250, 'LineStyle', 'none');
        colorbar;
        title(['Temperature Field at Time Step ', num2str(t)]);
        xlabel('x');
        ylabel('y');
    end
end

% Plot final temperature field
figure;
contourf(X, Y, T', 250, 'LineStyle', 'none');
colorbar;
title('Final Temperature Field');
xlabel('x');
ylabel('y');