function G = model_parameters(Nx, Ny, Nx_stim, Ny_stim, with_cells, Nm, lx, ly, Gg, dt, Tstop)
% G = model_parameters(Nx, Ny, Nx_stim, Ny_stim, with_cells, Nm, lx, ly, Gg, dt, Tstop)

% Set up discretization parameters
G.Nx = Nx;
G.Ny = Ny;
G.dt = dt; 
G.dt_ode = min(0.0005, dt); 
G.Tstop = Tstop;
G.nt = round(G.dt/G.dt_ode); 
G.Nt = round(G.Tstop/G.dt);
G.N = G.Nx*G.Ny;

% Set up stimulation information
G.stim_time = 1; % ms
G.stim_start_x = 1;
G.stim_start_y = 1;
G.stim_length_x = Nx_stim;
G.stim_length_y = Ny_stim;

% Specify cellular tissue part of the domain
G.with_myocyte = with_cells;
G.Nm = Nm;

% Define ionic model (human ventricular base model)
G.param_filename = 'membrane_model/parameterizations/C4_baseline.mat';
G.ion_model = 'Base_model';
G.Cm = 1; % uF/cm^2
[~, state_names] = base_model_init_states(G.param_filename);
G.V_idx = find(strcmp(state_names, 'V_m'));
G.init_states = @base_model_init_states;
G.init_parameters = @base_model_init_parameters;
[~, G.param_names] = G.init_parameters(G.param_filename);
G.rhs = @base_model_rhs_vectorized;


% Gap junction resistance
G.Rg_mx = (1/Gg)*ones(G.Nx-1, G.Ny); % kOhm
G.Rg_my = (1/Gg)*ones(G.Nx, G.Ny-1); % kOhm

% Adjust conductances at the cell collection boundary
with_idx = zeros(G.N, 1);
with_idx(G.with_myocyte) = 1;
with_idx = reshape(with_idx, G.Nx, G.Ny);
for i=1:G.Nx
    for j=1:G.Ny
        if with_idx(i,j) == 0
            if i > 1
                G.Rg_mx(i-1,j) = inf;
            end
            if i < G.Nx
                G.Rg_mx(i,j) = inf;
            end
            if j > 1
                G.Rg_my(i,j-1) = inf;
            end
            if j < G.Ny
                G.Rg_my(i,j) = inf;
            end
        end
    end
end

% Cell size
G.lx = lx; % cm
G.ly = ly; % cm
G.lz = ly; % cm
G.Am = zeros(G.Nx, G.Ny);
G.Am(G.with_myocyte) = 2*4*G.lx*G.ly;

% Volume fractions of cells and extracellular space
G.delta_m = zeros(G.Nx, G.Ny);
G.delta_m(G.with_myocyte) = 0.8;
G.delta_e = 1-G.delta_m;

% Conductivities and specific membrane capacitance
G.Cm = 1;        % uF/cm^2
G.sigma_e = 15;  % mS/cm
G.sigma_m = 8.2; % mS/cm


end
