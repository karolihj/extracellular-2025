% Run a Kirchhoff's network model simulation of a pancreatic islet of
% beta-cells
clear all
restoredefaultpath;
addpath('Riz2014/') % Membrane model

% Simulation options
p = 0.05; % Parameter variation fraction

% Set up simulation time
Tstop = 700e3; % Total simulation time (ms)

% Set up number of cells
N = 1000;

% Load model parameters
load('G_N1000.mat', 'G')
G.Tstop = Tstop;
G.Nt = round(G.Tstop/G.dt);
[G.Mm, G.Me] = set_up_conductances(G);

% Specify which parameters should be varied between cells and specify files 
% in which the parameter values for the individual cells are specified
if p > 0
    G.parameters_from_file = {'g_KATP_hat', 'g_CaL', 'g_Na'};
    G.parameter_files = {sprintf('cell_property_distributions/g_KATP_hat_N%d_p%d.txt', N, 100*p), ...
        sprintf('cell_property_distributions/g_CaL_N%d_p%d.txt', N, 100*p),  ...
        sprintf('cell_property_distributions/g_Na_N%d_p%d.txt', N, 100*p)};
end


% Run a simulation
[V, ~, U] = solve_system(G);

% Plot solution
t = 0:G.DT:G.Tstop;
U_max = mean(abs(U),2);
[~, sort_idx] = sort(U_max);
plot_idx = [sort_idx(end), sort_idx(1)];

subplot(2,1,1)
plot(t, V(plot_idx,:), 'linewidth', 2)
ylabel('v (mV)')

subplot(2,1,2)
plot(t, U(plot_idx,:)*1e3, 'linewidth', 2)
ylabel('u_e (\muV)')
xlabel('t (ms)')




