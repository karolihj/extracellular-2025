% Run a simulation of a simplified neuron modeled by the cable equation
clear all

% Specify simulation options
Na_move_percentage = 10;
use_stim = 0;

% Set up cable equation parameters
Tstop = 10000; % Total simulation time (in ms)
L = 1000e-4;   % Cell length (in cm)
d = 10e-4;     % Cell diameter (in cm)
sigma_i = 8.2; % Intracellular conductivity (in mS/cm)
sigma_e = 3;   % Extracellular conductivity (in mS/cm)
Cm = 1;        % Specific membrane capacitanceuF/cm^2
N = 1000;      % Number of discrete compartments
dx = L/N;      % Length of each compartment (in cm)
eta = d*sigma_i/4;

% Set up membrane model parameters
[param, param_names] = Masoli2015_init_parameters();
Np = length(param);
param = repmat(param, N, 1);

% Adjust g_Na
idx = round(N/2);
gNa_idx = find(strcmp(param_names, 'g_Na'));
gNa_bar = param(gNa_idx);
gNa_small = (1 - Na_move_percentage/100)*gNa_bar;
gNa_large = gNa_small + (Na_move_percentage/100)*(100/4)*gNa_bar;
param((0:N-1)*Np + gNa_idx) = gNa_small;
param((idx-20:idx+19)*Np + gNa_idx) = gNa_large;

% Set up stimulation
if use_stim
    stim_idx = find(strcmp(param_names, 'stim_amplitude'));
    stim_start_idx = find(strcmp(param_names, 'stim_start'));
    param((0:N-1)*Np + stim_start_idx) = 1;
    param((idx-1)*Np + stim_idx) = 3;
    param((idx)*Np + stim_idx) = 3;
    param((idx-2)*Np + stim_idx) = 3;
    param((idx+1)*Np + stim_idx) = 3;
end

% Set up initial conditions
[states, state_names] = Masoli2015_init_states();
Ns = length(states);
V_idx = find(strcmp(state_names, 'v'));
states = repmat(states, N, 1);

% Set up G matrix
G = spdiags([-1; -2*ones(N-2,1); -1], 0, N, N);
G = G + spdiags(ones(N, 1), 1, N, N);
G = G + spdiags(ones(N, 1), -1, N, N);
G = G*eta/(Cm*dx^2);

% Set up solver options
S_pattern = set_up_sparsity_pattern(N, Ns, V_idx);
if use_stim
    options = odeset('JPattern', S_pattern, 'MaxStep', 0.01);
else
    options = odeset('JPattern', S_pattern);
end

% Perform a cable model simulation
P.parameters = param; P.N = N; P.G = G;
[T, S] = ode15s(@system_rhs, [0, Tstop], states, options, P);
V = S(:, V_idx:Ns:end);

% Extract the current sources for comuputing the extracellular potential
Im = G*V';        % Transmenbrane current density
c = (dx*d*pi)*Im; % Current source

% Set up spatial points for computing the extracellular potential
x = (dx/2:dx:L-dx/2);
z = zeros(N, 1);
y_cell = zeros(N, 1);
y_e = 6e-4*ones(N, 1);

% Compute the extracellular potential
ue = zeros(N, length(T));
for n=1:N
    for k=1:N
        ue(n,:) = ue(n,:) + (1/(4*pi*sigma_e))*c(k,:)/(norm([x(n);y_e(n);z(n)]-[x(k);y_cell(k);z(k)]));
    end
end


% Plot the membrane potential and the extracellular potential in two points
plot_idx = [idx + 1, idx + round(N/6)];
subplot(2,1,1)
plot(T, V(:, plot_idx), 'LineWidth', 2);
ylabel('v (mV)')

subplot(2,1,2)
plot(T, ue(plot_idx,:), 'LineWidth', 2);
ylabel('u_e (\muV)')
xlabel('t (ms)')




