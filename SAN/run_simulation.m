% Run a bidomain model simulation of a sinoatrial node tissue sample
clear all
addpath('Fabbri2017/') % Load membrane model

% Simulation options
p = 0.05; % Parameter variation fraction

% Set up simulation time
Tstop = 30000; % Total simulation time (in ms)

% Set up spatial discretization
dx = 0.002; % cm
dy = 0.002; % cm
dz = 0.002; % cm

% Set up domain sizes
Lx = 0.05; % cm
Ly = 0.02; % cm
Lz = 0.01; % cm
Lb = 0.01; % Extracellular bath length (in cm)

% Generate model parameters
G = model_parameters(Lx, Ly, Lz, dx, dy, dz, Lb, Tstop);
[G.Mm, G.Me] = set_up_conductivities(G);
G.vary_percentage = p;

% Run a bidomain model simulation
[V, states, U] = solve_system(G);

% Plot the solution
t = 0:G.DT:G.Tstop;
plot_idx = [((G.Nz-round(G.Lb/G.dz)-1)*G.Ny + round((G.Lb+G.Ly*0.75)/G.dy))*G.Nx + round((G.Lb+0.25*G.Lx)/G.dx), ...
    ((G.Nz-round(G.Lb/G.dz)-1)*G.Ny + round((G.Lb+G.Ly*0.25)/G.dy))*G.Nx + round((G.Lb+0.75*G.Lx)/G.dx)];

subplot(2,1,1)
plot(t, V(plot_idx,:), 'linewidth', 2)
ylabel('v (mV)')

subplot(2,1,2)
plot(t, U(plot_idx,:)*1e3, 'linewidth', 2)
ylabel('u_e (\muV)')
xlabel('t (ms)')





