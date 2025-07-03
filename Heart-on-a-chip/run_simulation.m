% Run a Kirchhoff's network model simulation of a heart-on-a-chip with
% hiPSC-CMs
restoredefaultpath;
addpath('membrane_model/')

% Simulation options
stimulate = 1;
vary_percentage = 0.1; % Parameter variation fraction

% Set up temporal parameters
Tstop = 12000; % Total simulation time (ms)
dt = 0.1;      % Time step (ms)

% Load cell locations
load('locations/C4_cells.mat', 'lx', 'ly', 'Nx', 'Ny', 'N', 'Nm', 'with_cells')

% Set up stimulation location
Nx_stim = round(350/lx);
Ny_stim = Ny;

% Set up additional model parameters
Gg = 180e-6; % Gap junction conductance (in mS)
G = model_parameters(Nx, Ny, Nx_stim, Ny_stim, with_cells, Nm, lx*1e-4, ly*1e-4, Gg, dt, Tstop);
[G.Mm, G.Me] = set_up_conductances(G);

if stimulate
    G.stim_amp = 15;
else
    G.stim_amp = 0;
end
G.vary_percentage = vary_percentage;

% Set up mesh
mesh = set_up_mesh(G);

% Run a simulation
[V, U] = solve_system(G, mesh);

% Plot solution
t = 0:G.dt:G.Tstop;
plot_idx = G.with_myocyte(858);

subplot(2,1,1)
plot(t, V(plot_idx,:), 'linewidth', 2)
ylabel('v (mV)')

subplot(2,1,2)
plot(t, U(plot_idx,:), 'linewidth', 2)
ylabel('u_e (mV)')
xlabel('t (ms)')