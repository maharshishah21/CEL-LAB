% Define heat plate dimensions and properties
Nx = 25; % Number of grid points in x-direction
Ny = 25; % Number of grid points in y-direction
Lx = 1; % Length of the plate in x-direction (m)
Ly = 1; % Length of the plate in y-direction (m)
Lt = 750; % Number of time steps
dx = Lx / (Nx - 1); % Grid spacing in x-direction
dy = Ly / (Ny - 1); % Grid spacing in y-direction
dt = 0.1; % Time step size
alpha = 0.0001; % Thermal diffusivity (m^2/s)
T_high = 50; % Elevated temperature
T_low = 10; % Lower temperature

% Set initial temperature (half plate at elevated temperature, half at lower temperature)
T = T_low * ones(Ny, Nx);
T(round(Ny/2):end, :) = T_high;

% Create a figure with enough subplots for all time steps
figure;

% Plot initial temperature field 
subplot(2, 2, 1);
[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));
contourf(X, Y, T', 250, 'LineStyle', 'none');
colorbar;
title('Initial Temperature Field');

% Main simulation loop
for t = 1:Lt
    T_old = T;
    for i = 2:Nx-1
        for j = 2:Ny-1
            T(i, j) = T_old(i, j) + alpha * dt * ((T_old(i+1, j) - 2*T_old(i, j) + T_old(i-1, j)) / dx^2 + ... 
                (T_old(i, j+1) - 2*T_old(i, j) + T_old(i, j-1)) / dy^2);        
        end
    end
end

% Plot final temperature field
subplot(2, 2, 2);
contourf(X, Y, T', 250, 'LineStyle', 'none');
colorbar;
title('Final Temperature Field');