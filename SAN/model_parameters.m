function G = model_parameters(Lx, Ly, Lz, dx, dy, dz, Lb, Tstop)
% G = model_parameters(Lx, Ly, Lz, dx, dy, dz, Lb, Tstop)

% Set up geometry
G.Lx = Lx; G.Ly = Ly; G.Lz = Lz; G.Lb = Lb; G.dx = dx; G.dy = dy; G.dz = dz;
G.Nx = round((2*Lb+Lx)/dx) + 1;
G.Ny = round((2*Lb+Ly)/dy) + 1;
G.Nz = round((2*Lb+Lz)/dz) + 1;
G.N = G.Nx*G.Ny*G.Nz;
G.include_ue = 1;

% Remove cells from the extracellular bath
G.with_cell = zeros(G.N, 1);
counter = 1;
for k=1:G.Nz
    for j=1:G.Ny
        for i=1:G.Nx
            if (G.dx*(i-1) >= G.Lb) && (G.dx*(i-1) <= G.Lx+G.Lb) ...
                    && (G.dy*(j-1) >= G.Lb) && (G.dy*(j-1) <= G.Ly+G.Lb) ...
                    && (G.dz*(k-1) >= G.Lb) && (G.dz*(k-1) <= G.Lz+G.Lb)
                G.with_cell(counter) = 1;
            end
            counter = counter + 1;
        end
    end
end
G.with_cell = find(G.with_cell);
G.Nm = length(G.with_cell);


% Set up temporal discretization parameters
G.dt = 0.1;                   % Time step for the pde step 
G.dt_ode = min(0.01, G.dt);   % Time step for the ode step 
G.Tstop = Tstop;              % Total simulation time
G.nt = round(G.dt/G.dt_ode);  % Number of ode steps per pde step 
G.Nt = round(G.Tstop/G.dt);   % Number of pde steps 
G.DT = max(1, G.dt);          % Time step for saving the solution

% Define membrane model
G.init_states = @Fabbri2017_init_states;
G.init_parameters = @Fabbri2017_init_parameters;
[~, state_names] = G.init_states();
G.V_idx = find(strcmp(state_names, 'V_ode'));
G.rhs = @Fabbri2017_rhs_vectorized;

% Gap junction resistance
G.Rg_mx = 1e6*ones(G.Nx-1, G.Ny, G.Nz); % kOhm
G.Rg_my = 1e6*ones(G.Nx, G.Ny-1, G.Nz); % kOhm
G.Rg_mz = 1e6*ones(G.Nx, G.Ny, G.Nz-1); % kOhm

% Remove gap junction resistances at cellular tissue boundary
with_idx = zeros(G.N, 1);
with_idx(G.with_cell) = 1;
with_idx = reshape(with_idx, G.Nx, G.Ny, G.Nz);
for i=1:G.Nx
    for j=1:G.Ny
        for k=1:G.Nz
            if with_idx(i,j,k) == 0
                if i > 1
                    G.Rg_mx(i-1,j,k) = inf;
                end
                if i < G.Nx
                    G.Rg_mx(i,j,k) = inf;
                end
                if j > 1
                    G.Rg_my(i,j-1,k) = inf;
                end
                if j < G.Ny
                    G.Rg_my(i,j,k) = inf;
                end
                if k > 1
                    G.Rg_mz(i,j,k-1) = inf;
                end
                if k < G.Nz
                    G.Rg_mz(i,j,k) = inf;
                end
            end
        end
    end
end

% Cell size
G.lx = 67e-4;  % cm
G.ly = 7.8e-4; % cm
G.lz = 7.8e-4; % cm
G.chi = zeros(G.N, 1);
G.chi(G.with_cell) = (57/1e6)/(G.lx*G.ly*G.lz);

% Volume fractions of cells and extracellular space
G.delta_m = zeros(G.N, 1);
G.delta_m(G.with_cell) = 0.8;
G.delta_m = reshape(G.delta_m, G.Nx, G.Ny, G.Nz);
G.delta_e = 1-G.delta_m;

% Conductivities and capacitance
G.Cm = 1;        % uF/cm^2
G.sigma_e = 15;  % mS/cm
G.sigma_m = 8.2; % mS/cm

end
